import { mysqlTable, int, varchar, text, timestamp, index } from 'drizzle-orm/mysql-core';

export const dockingQueue = mysqlTable('docking_queue', {
  id: int('id').primaryKey().autoincrement(),
  analogId: int('analog_id').notNull(),
  target: varchar('target', { length: 100 }).notNull(), // 'NMDA', '5-HT2A', 'D2'
  status: varchar('status', { length: 50 }).notNull().default('pending'), // 'pending', 'running', 'completed', 'failed'
  priority: int('priority').notNull().default(5), // 1-10, higher = more important
  
  // Results
  bindingAffinity: varchar('binding_affinity', { length: 50 }),
  dockingScore: int('docking_score'),
  ligandPDB: text('ligand_pdb'),
  receptorPDB: text('receptor_pdb'),
  
  // Metadata
  errorMessage: text('error_message'),
  startedAt: timestamp('started_at'),
  completedAt: timestamp('completed_at'),
  createdAt: timestamp('created_at').notNull().defaultNow(),
  updatedAt: timestamp('updated_at').notNull().defaultNow().onUpdateNow(),
}, (table) => ({
  analogIdIdx: index('analog_id_idx').on(table.analogId),
  statusIdx: index('status_idx').on(table.status),
  targetIdx: index('target_idx').on(table.target),
}));

export type DockingQueueEntry = typeof dockingQueue.$inferSelect;
export type NewDockingQueueEntry = typeof dockingQueue.$inferInsert;
