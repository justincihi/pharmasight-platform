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
  const venvPython = path.join(__dirname, '..', 'python_modules', 'venv', 'bin', 'python');
  const script = path.join(__dirname, '..', '..', 'python_services', 'analysis_microservice.py');
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
  child.on('exit', (code) => console.warn(`[MicroService] exited with code ${code}`));
  console.log('[MicroService] Analysis microservice starting on port 5000...');
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
    // Initialize autonomous research scheduler (async: loads cron from DB)
    initializeScheduler().catch(console.error);
  });
}

startServer().catch(console.error);
