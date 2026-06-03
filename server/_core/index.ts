import "dotenv/config";
import express from "express";
import { createServer } from "http";
import net from "net";
import { spawn } from "child_process";
import path from "path";
import { fileURLToPath } from "url";
import { createExpressMiddleware } from "@trpc/server/adapters/express";
import { registerOAuthRoutes } from "./oauth";
import { appRouter } from "../routers";
import { createContext } from "./context";
import { serveStatic, setupVite } from "./vite";
import { initializeScheduler } from "../initScheduler";
import { registerPlatformAPI } from "../platformAPI";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

/** Start the Flask analysis microservice in the background (port 5000) */
function startAnalysisMicroservice(): void {
  // In production (Cloud Run) Python/venv are not available — skip silently.
  // All Python-dependent procedures fall back to LLM or mock responses.
  if (process.env.NODE_ENV === 'production') {
    console.log('[MicroService] Production environment — skipping local Python microservice.');
    return;
  }
  // Check if port 5000 is already occupied before spawning
  const probe = net.createServer();
  probe.once('error', () => {
    // Port already in use — microservice is already running, skip
    console.log('[MicroService] Port 5000 already in use, skipping start.');
  });
  probe.once('listening', () => {
    probe.close(() => {
      // Port is free — start the microservice
      const venvPython = path.join(__dirname, '..', 'python_modules', 'venv', 'bin', 'python');
      const script = path.join(__dirname, '..', '..', 'python_services', 'analysis_microservice.py');
      try {
        const child = spawn(venvPython, [script], {
          env: { ...process.env, PORT: '5000' },
          stdio: ['ignore', 'pipe', 'pipe'],
          detached: false,
        });
        child.stdout?.on('data', (d: Buffer) => {
          const line = d.toString().trim();
          if (line) console.log('[MicroService]', line);
        });
        child.stderr?.on('data', (d: Buffer) => {
          const line = d.toString().trim();
          if (line && !line.includes('WARNING') && !line.includes('Serving Flask') && !line.includes('Press CTRL')) {
            console.error('[MicroService]', line);
          }
        });
        child.on('error', (err) => console.warn('[MicroService] Spawn error (non-fatal):', err.message));
        child.on('exit', (code) => console.warn(`[MicroService] exited with code ${code}`));
        console.log('[MicroService] Analysis microservice starting on port 5000...');
      } catch (err) {
        console.warn('[MicroService] Could not start Python microservice (non-fatal):', (err as Error).message);
      }
    });
  });
  probe.on('error', (err) => console.warn('[MicroService] Port probe error (non-fatal):', err.message));
  probe.listen(5000, '127.0.0.1');
}

function isPortAvailable(port: number): Promise<boolean> {
  return new Promise(resolve => {
    const server = net.createServer();
    server.listen(port, () => {
      server.close(() => resolve(true));
    });
    server.on("error", () => resolve(false));
  });
}

async function findAvailablePort(startPort: number = 3000): Promise<number> {
  for (let port = startPort; port < startPort + 20; port++) {
    if (await isPortAvailable(port)) {
      return port;
    }
  }
  throw new Error(`No available port found starting from ${startPort}`);
}

/** Mark any research runs left in 'running' state (from a crashed/restarted server) as failed */
async function cleanupStaleRuns(): Promise<void> {
  try {
    const { getDb } = await import('../db');
    const { researchRuns } = await import('../../drizzle/schema');
    const { eq, lt } = await import('drizzle-orm');
    const db = await getDb();
    if (!db) return;
    // Any run still 'running' after server boot is orphaned — mark as failed
    const staleThreshold = new Date(Date.now() - 30 * 60 * 1000); // 30 min ago
    const result = await db
      .update(researchRuns)
      .set({
        status: 'failed',
        errorMessage: 'Run interrupted: server restarted before completion.',
        completedAt: new Date(),
      })
      .where(eq(researchRuns.status, 'running'));
    const affected = (result as any).rowsAffected ?? (result as any)[0]?.affectedRows ?? 0;
    if (affected > 0) {
      console.log(`[Init] Cleaned up ${affected} stale research run(s) left in 'running' state.`);
    }
  } catch (e) {
    console.error('[Init] Failed to clean up stale runs:', e);
  }
}

async function startServer() {
  const app = express();
  const server = createServer(app);
  // Configure body parser with larger size limit for file uploads
  app.use(express.json({ limit: "50mb" }));
  app.use(express.urlencoded({ limit: "50mb", extended: true }));
  // OAuth callback under /api/oauth/callback
  registerOAuthRoutes(app);
  // Platform API for PharmaSight integration
  registerPlatformAPI(app);
  // tRPC API
  app.use(
    "/api/trpc",
    createExpressMiddleware({
      router: appRouter,
      createContext,
    })
  );
  // development mode uses Vite, production mode uses static files
  if (process.env.NODE_ENV === "development") {
    await setupVite(app, server);
  } else {
    serveStatic(app);
  }

  const preferredPort = parseInt(process.env.PORT || "3000");
  const port = await findAvailablePort(preferredPort);

  if (port !== preferredPort) {
    console.log(`Port ${preferredPort} is busy, using port ${port} instead`);
  }

  server.listen(port, () => {
    console.log(`Server running on http://localhost:${port}/`);
    // Start Python analysis microservice (metabolites, ADMET, docking)
    startAnalysisMicroservice();
    // Clean up any runs that were left in 'running' state from a previous server instance
    cleanupStaleRuns().catch(console.error);
    // Initialize autonomous research scheduler (async: loads cron from DB)
    initializeScheduler().catch(console.error);
  });
}

startServer().catch(console.error);
