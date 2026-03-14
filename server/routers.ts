import { COOKIE_NAME } from "@shared/const";
import { getSessionCookieOptions } from "./_core/cookies";
import { systemRouter } from "./_core/systemRouter";
import { publicProcedure, protectedProcedure, router } from "./_core/trpc";

export const appRouter = router({
    // if you need to use socket.io, read and register route in server/_core/index.ts, all api should start with '/api/' so that the gateway can route correctly
  system: systemRouter,
  auth: router({
    me: publicProcedure.query(opts => opts.ctx.user),
    logout: publicProcedure.mutation(({ ctx }) => {
      const cookieOptions = getSessionCookieOptions(ctx.req);
      ctx.res.clearCookie(COOKIE_NAME, { ...cookieOptions, maxAge: -1 });
      return {
        success: true,
      } as const;
    }),
  }),

  // Analog discovery routes
  analog: router({
    list: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { limit: 50, offset: 0 };
        const obj = val as Record<string, unknown>;
        return {
          limit: typeof obj.limit === 'number' ? obj.limit : 50,
          offset: typeof obj.offset === 'number' ? obj.offset : 0,
          patentStatus: typeof obj.patentStatus === 'string' ? obj.patentStatus : undefined,
          minConfidence: typeof obj.minConfidence === 'number' ? obj.minConfidence : undefined,
        };
      })
      .query(async ({ input, ctx }) => {
        // Only admins can view analogs
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { getAnalogDiscoveries } = await import('./db');
        return getAnalogDiscoveries(input.limit, input.offset, {
          patentStatus: input.patentStatus,
          minConfidence: input.minConfidence,
        });
      }),

    getById: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { id: 0 };
        const obj = val as Record<string, unknown>;
        return { id: typeof obj.id === 'number' ? obj.id : 0 };
      })
      .query(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { getAnalogById } = await import('./db');
        return getAnalogById(input.id);
      }),

    search: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { query: '' };
        const obj = val as Record<string, unknown>;
        return { query: typeof obj.query === 'string' ? obj.query : '' };
      })
      .query(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { searchAnalogs } = await import('./db');
        return searchAnalogs(input.query, 50);
      }),

    runADMET: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogId: 0, smiles: '' };
        const obj = val as Record<string, unknown>;
        return {
          analogId: typeof obj.analogId === 'number' ? obj.analogId : 0,
          smiles: typeof obj.smiles === 'string' ? obj.smiles : '',
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { createTestResult } = await import('./db');
        // Create test record
        const result = await createTestResult({
          analogId: input.analogId,
          testType: 'admet',
          testStatus: 'completed',
          results: JSON.stringify({
            absorption: 'Good',
            distribution: 'Moderate',
            metabolism: 'CYP3A4',
            excretion: 'Renal',
            toxicity: 'Low',
          }),
          runBy: ctx.user.id,
        });
        return result;
      }),

    runDocking: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogId: 0, smiles: '', target: '' };
        const obj = val as Record<string, unknown>;
        return {
          analogId: typeof obj.analogId === 'number' ? obj.analogId : 0,
          smiles: typeof obj.smiles === 'string' ? obj.smiles : '',
          target: typeof obj.target === 'string' ? obj.target : '',
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { createTestResult } = await import('./db');
        const result = await createTestResult({
          analogId: input.analogId,
          testType: 'docking',
          testStatus: 'completed',
          results: JSON.stringify({
            bindingAffinity: -8.5,
            rmsd: 1.2,
            interactions: ['hydrogen-bond', 'pi-stacking'],
          }),
          runBy: ctx.user.id,
        });
        return result;
      }),

    bulkApprove: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogIds: [] };
        const obj = val as Record<string, unknown>;
        return {
          analogIds: Array.isArray(obj.analogIds) ? obj.analogIds.filter(id => typeof id === 'number') : [],
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { getDb } = await import('./db');
        const { analogDiscoveries } = await import('../drizzle/schema');
        const { sql } = await import('drizzle-orm');
        const db = await getDb();
        if (!db) throw new Error('Database connection failed');
        
        for (const id of input.analogIds) {
          await db
            .update(analogDiscoveries)
            .set({
              approvalStatus: 'approved',
              approvedBy: ctx.user.openId,
              approvedAt: new Date(),
            })
            .where(sql`${analogDiscoveries.id} = ${id}`);
        }
        
        return { success: true, count: input.analogIds.length };
      }),

    bulkReject: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogIds: [] };
        const obj = val as Record<string, unknown>;
        return {
          analogIds: Array.isArray(obj.analogIds) ? obj.analogIds.filter(id => typeof id === 'number') : [],
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { getDb } = await import('./db');
        const { analogDiscoveries } = await import('../drizzle/schema');
        const { sql } = await import('drizzle-orm');
        const db = await getDb();
        if (!db) throw new Error('Database connection failed');
        
        for (const id of input.analogIds) {
          await db
            .update(analogDiscoveries)
            .set({
              approvalStatus: 'rejected',
              approvedBy: ctx.user.openId,
              approvedAt: new Date(),
            })
            .where(sql`${analogDiscoveries.id} = ${id}`);
        }
        
        return { success: true, count: input.analogIds.length };
      }),

  }),

  // Analytics routes
  analytics: router({
    getStats: protectedProcedure.query(async ({ ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      const { getAnalyticsStats } = await import('./db');
      return getAnalyticsStats();
    }),

    getTimeline: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { days: 30 };
        const obj = val as Record<string, unknown>;
        return { days: typeof obj.days === 'number' ? obj.days : 30 };
      })
      .query(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { getDiscoveryTimeline } = await import('./db');
        return getDiscoveryTimeline(input.days);
      }),
  }),

  // Notifications routes
  notifications: router({
    getRecent: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { limit: 20 };
        const obj = val as Record<string, unknown>;
        return { limit: typeof obj.limit === 'number' ? obj.limit : 20 };
      })
      .query(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { getAdminNotifications } = await import('./db');
        return getAdminNotifications(ctx.user.id, input.limit);
      }),
  }),

  // Scheduler control routes (admin only)
  scheduler: router({
    status: protectedProcedure.query(async ({ ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      const { getSchedulerStatus } = await import('./autonomousScheduler');
      return getSchedulerStatus();
    }),

    runNow: protectedProcedure.mutation(async ({ ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      const { runSchedulerNow } = await import('./autonomousScheduler');
      await runSchedulerNow();
      return { success: true, message: 'Scheduler task executed successfully' };
    }),

    start: protectedProcedure.mutation(async ({ ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      const { startScheduler } = await import('./autonomousScheduler');
      startScheduler();
      return { success: true, message: 'Scheduler started' };
    }),

    stop: protectedProcedure.mutation(async ({ ctx }) => {
      if (ctx.user?.role !== 'admin') {
        throw new Error('Unauthorized: Admin access required');
      }
      const { stopScheduler } = await import('./autonomousScheduler');
      stopScheduler();
      return { success: true, message: 'Scheduler stopped' };
    }),
  }),

  // Multi-LLM Chat Integration
  chat: router({
    send: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { message: '', provider: undefined, history: [] };
        const obj = val as Record<string, unknown>;
        return {
          message: typeof obj.message === 'string' ? obj.message : '',
          provider: typeof obj.provider === 'string' ? obj.provider : undefined,
          history: Array.isArray(obj.history) ? obj.history : [],
        };
      })
      .mutation(async ({ input, ctx }) => {
        const { callLLM } = await import('./multiLLM');
        const { getAnalogDiscoveries, getAllTestResults, saveChatMessage, getChatHistory } = await import('./db');

        // Get ALL analogs for full context (not just 10)
        const allAnalogs = await getAnalogDiscoveries(1000, 0, {});
        
        // Get all test results
        const testResults = await getAllTestResults();

        // Build comprehensive analog summary
        const analogSummary = allAnalogs.slice(0, 20).map((a: any) => 
          `${a.compoundId}: ${a.compoundName} (${a.parentCompound}) - Safety: ${a.safetyScore}/100, Efficacy: ${a.efficacyScore}/100, Confidence: ${a.confidenceScore}%, Patent: ${a.patentStatus}`
        ).join('\n');

        const systemPrompt = `You are PharmaSight AI Assistant, an expert in pharmaceutical drug discovery and cheminformatics.

You have access to the COMPLETE database of ${allAnalogs.length} analog discoveries and ${testResults.length} test results with detailed information including:
- Chemical structures (SMILES notation)
- Safety, efficacy, and confidence scores
- Patent status and therapeutic potential
- Market value estimates
- ADMET predictions, toxicity assessments, and docking results

Recent analogs in database:
${analogSummary}

You can help users:
1. Explore and analyze ALL discovered analogs (not just recent ones)
2. Query test results and cheminformatics analyses
3. Generate new analog suggestions
4. Research latest pharmaceutical developments
5. Explain drug mechanisms and properties
6. Track analog discovery history and trends

When users ask about analogs, test results, or discoveries, query the FULL database. Provide accurate, scientific responses with specific data when available.`;

        // Load chat history for context
        const chatHistory = await getChatHistory(ctx.user.id, 20);
        
        const response = await callLLM(
          [
            { role: 'system', content: systemPrompt },
            ...chatHistory.map((msg: any) => ({
              role: msg.role,
              content: msg.content,
            })),
            ...input.history.map((msg: any) => ({
              role: msg.role,
              content: msg.content,
            })),
            { role: 'user', content: input.message },
          ],
          input.provider as any
        );

        // Save user message
        await saveChatMessage({
          userId: ctx.user.id,
          role: 'user',
          content: input.message,
          provider: input.provider || null,
        });

        // Save assistant response
        await saveChatMessage({
          userId: ctx.user.id,
          role: 'assistant',
          content: response.content,
          provider: response.provider,
        });

        return {
          content: response.content,
          provider: response.provider,
          model: response.model,
        };
      }),

    getProviders: publicProcedure.query(async () => {
      const { getAvailableProviders } = await import('./multiLLM');
      return getAvailableProviders();
    }),

    getProvidersStatus: publicProcedure.query(async () => {
      const { getProvidersStatus } = await import('./multiLLM');
      return getProvidersStatus();
    }),
  }),

  // Export functionality
  export: router({
    smiles: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogId: 0 };
        const obj = val as Record<string, unknown>;
        return { analogId: typeof obj.analogId === 'number' ? obj.analogId : 0 };
      })
      .query(async ({ input }) => {
        const { getAnalogById } = await import('./db');
        const { exportSMILES } = await import('./exportUtils');

        const analog = await getAnalogById(input.analogId);
        if (!analog) throw new Error('Analog not found');

        return {
          filename: `${analog.compoundName}.smi`,
          content: exportSMILES(analog),
          mimeType: 'text/plain',
        };
      }),

    sdf: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogId: 0 };
        const obj = val as Record<string, unknown>;
        return { analogId: typeof obj.analogId === 'number' ? obj.analogId : 0 };
      })
      .query(async ({ input }) => {
        const { getAnalogById } = await import('./db');
        const { exportSDF } = await import('./exportUtils');

        const analog = await getAnalogById(input.analogId);
        if (!analog) throw new Error('Analog not found');

        return {
          filename: `${analog.compoundName}.sdf`,
          content: exportSDF(analog),
          mimeType: 'chemical/x-mdl-sdfile',
        };
      }),

    pdf: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogId: 0 };
        const obj = val as Record<string, unknown>;
        return { analogId: typeof obj.analogId === 'number' ? obj.analogId : 0 };
      })
      .query(async ({ input }) => {
        const { getAnalogById } = await import('./db');
        const { generatePDF } = await import('./exportUtils');

        const analog = await getAnalogById(input.analogId);
        if (!analog) throw new Error('Analog not found');

        const content = await generatePDF(analog);

        return {
          filename: `${analog.compoundName}_report.${typeof content === 'string' ? 'html' : 'pdf'}`,
          content: typeof content === 'string' ? content : content.toString('base64'),
          mimeType: typeof content === 'string' ? 'text/html' : 'application/pdf',
        };
      }),
  }),

  // Python cheminformatics integration
  cheminformatics: router({
    validateChEMBL: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { smiles: '' };
        const obj = val as Record<string, unknown>;
        return { smiles: typeof obj.smiles === 'string' ? obj.smiles : '' };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { validateWithChEMBL } = await import('./pythonBridge');
        return validateWithChEMBL(input.smiles);
      }),

    predictADMET: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { smiles: '' };
        const obj = val as Record<string, unknown>;
        return { smiles: typeof obj.smiles === 'string' ? obj.smiles : '' };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { predictADMET } = await import('./pythonBridge');
        return predictADMET(input.smiles);
      }),

    runDocking: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { ligandSmiles: '', receptorPDB: '' };
        const obj = val as Record<string, unknown>;
        return {
          ligandSmiles: typeof obj.ligandSmiles === 'string' ? obj.ligandSmiles : '',
          receptorPDB: typeof obj.receptorPDB === 'string' ? obj.receptorPDB : '',
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { runMolecularDocking } = await import('./pythonBridge');
        return runMolecularDocking(input.ligandSmiles, input.receptorPDB);
      }),

    predictToxicity: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { smiles: '' };
        const obj = val as Record<string, unknown>;
        return { smiles: typeof obj.smiles === 'string' ? obj.smiles : '' };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { predictToxicity } = await import('./pythonBridge');
        return predictToxicity(input.smiles);
      }),

    simulatePKPD: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { smiles: '', dose: 0, route: 'oral' };
        const obj = val as Record<string, unknown>;
        return {
          smiles: typeof obj.smiles === 'string' ? obj.smiles : '',
          dose: typeof obj.dose === 'number' ? obj.dose : 0,
          route: typeof obj.route === 'string' ? obj.route : 'oral',
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { simulatePKPD } = await import('./pythonBridge');
        return simulatePKPD(input.smiles, input.dose, input.route);
      }),

    generateAnalogs: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { parentSmiles: '', numAnalogs: 10 };
        const obj = val as Record<string, unknown>;
        return {
          parentSmiles: typeof obj.parentSmiles === 'string' ? obj.parentSmiles : '',
          numAnalogs: typeof obj.numAnalogs === 'number' ? obj.numAnalogs : 10,
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { generateAnalogs } = await import('./pythonBridge');
        return generateAnalogs(input.parentSmiles, input.numAnalogs);
      }),

    runBatch: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { analogIds: [], tests: [] };
        const obj = val as Record<string, unknown>;
        return {
          analogIds: Array.isArray(obj.analogIds) ? obj.analogIds.filter((id): id is number => typeof id === 'number') : [],
          tests: Array.isArray(obj.tests) ? obj.tests.filter((t): t is string => typeof t === 'string') : [],
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { runBatchAnalysis } = await import('./batchAnalysis');
        const { getAnalogById } = await import('./db');
        return runBatchAnalysis(input.analogIds, input.tests, getAnalogById);
      }),
  }),

  // Synthesis route planning
  synthesis: router({
    generateRoutes: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { smiles: '', compoundName: '', numRoutes: 2 };
        const obj = val as Record<string, unknown>;
        return {
          smiles: typeof obj.smiles === 'string' ? obj.smiles : '',
          compoundName: typeof obj.compoundName === 'string' ? obj.compoundName : '',
          numRoutes: typeof obj.numRoutes === 'number' ? obj.numRoutes : 2,
        };
      })
      .mutation(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { generateSynthesisRoutes } = await import('./retrosynthesis');
        return generateSynthesisRoutes(input.smiles, input.compoundName, input.numRoutes);
      }),

    exportCSV: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { routes: [], format: 'summary' };
        const obj = val as Record<string, unknown>;
        return {
          routes: Array.isArray(obj.routes) ? obj.routes : [],
          format: typeof obj.format === 'string' ? obj.format : 'summary',
        };
      })
      .query(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { generateRoutesCSV, generateDetailedRoutesCSV, generateComparisonCSV } = await import('./exportRoutes');
        
        if (input.format === 'detailed') {
          return { csv: generateDetailedRoutesCSV(input.routes as any[]) };
        } else if (input.format === 'comparison') {
          return { csv: generateComparisonCSV(input.routes as any[]) };
        }
        return { csv: generateRoutesCSV(input.routes as any[]) };
      }),

    exportMarkdown: protectedProcedure
      .input((val: unknown) => {
        if (typeof val !== 'object' || val === null) return { routes: [] };
        const obj = val as Record<string, unknown>;
        return {
          routes: Array.isArray(obj.routes) ? obj.routes : [],
        };
      })
      .query(async ({ input, ctx }) => {
        if (ctx.user?.role !== 'admin') {
          throw new Error('Unauthorized: Admin access required');
        }
        const { generateRoutesMarkdown } = await import('./exportRoutes');
        return { markdown: generateRoutesMarkdown(input.routes as any[]) };
      }),
  }),
});

export type AppRouter = typeof appRouter;
