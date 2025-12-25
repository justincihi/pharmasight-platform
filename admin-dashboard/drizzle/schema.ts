import { int, mysqlEnum, mysqlTable, text, timestamp, varchar } from "drizzle-orm/mysql-core";

/**
 * Core user table backing auth flow.
 * Extend this file with additional tables as your product grows.
 * Columns use camelCase to match both database fields and generated types.
 */
export const users = mysqlTable("users", {
  /**
   * Surrogate primary key. Auto-incremented numeric value managed by the database.
   * Use this for relations between tables.
   */
  id: int("id").autoincrement().primaryKey(),
  /** Manus OAuth identifier (openId) returned from the OAuth callback. Unique per user. */
  openId: varchar("openId", { length: 64 }).notNull().unique(),
  name: text("name"),
  email: varchar("email", { length: 320 }),
  loginMethod: varchar("loginMethod", { length: 64 }),
  role: mysqlEnum("role", ["user", "admin"]).default("user").notNull(),
  createdAt: timestamp("createdAt").defaultNow().notNull(),
  updatedAt: timestamp("updatedAt").defaultNow().onUpdateNow().notNull(),
  lastSignedIn: timestamp("lastSignedIn").defaultNow().notNull(),
});

export type User = typeof users.$inferSelect;
export type InsertUser = typeof users.$inferInsert;

/**
 * Analog discoveries table - stores all discovered pharmaceutical analogs
 */
export const analogDiscoveries = mysqlTable("analog_discoveries", {
  id: int("id").autoincrement().primaryKey(),
  compoundId: varchar("compound_id", { length: 128 }).notNull().unique(),
  compoundName: varchar("compound_name", { length: 255 }).notNull(),
  parentCompound: varchar("parent_compound", { length: 255 }).notNull(),
  smiles: text("smiles").notNull(),
  
  // Scores and metrics
  confidenceScore: int("confidence_score").notNull(), // 0-100
  similarityScore: int("similarity_score").notNull(), // 0-100
  safetyScore: int("safety_score").notNull(), // 0-100
  efficacyScore: int("efficacy_score").notNull(), // 0-100
  drugLikenessScore: int("drug_likeness_score").notNull(), // 0-100
  
  // Patent and regulatory
  patentStatus: mysqlEnum("patent_status", ["patent-free", "patent-opportunity", "patented", "unknown"]).notNull(),
  patentNumbers: text("patent_numbers"), // JSON array of patent numbers
  fdaStatus: varchar("fda_status", { length: 64 }), // approved, investigational, etc.
  
  // Market and therapeutic
  marketValue: varchar("market_value", { length: 64 }), // e.g., "$45M"
  therapeuticPotential: text("therapeutic_potential"),
  keyDifferences: text("key_differences"),
  mechanismOfAction: text("mechanism_of_action"),
  
  // Chemical properties
  molecularWeight: varchar("molecular_weight", { length: 64 }),
  logP: varchar("log_p", { length: 64 }),
  hBondDonors: int("h_bond_donors"),
  hBondAcceptors: int("h_bond_acceptors"),
  
  // External database IDs
  pubchemCid: varchar("pubchem_cid", { length: 64 }),
  chemblId: varchar("chembl_id", { length: 64 }),
  
  // Discovery metadata
  discoveredBy: varchar("discovered_by", { length: 128 }).notNull(), // "autonomous-engine" or user ID
  discoveryMethod: varchar("discovery_method", { length: 128 }),
  discoveredAt: timestamp("discovered_at").defaultNow().notNull(),
  
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type AnalogDiscovery = typeof analogDiscoveries.$inferSelect;
export type InsertAnalogDiscovery = typeof analogDiscoveries.$inferInsert;

/**
 * Cheminformatics test results table - stores results from various analyses
 */
export const testResults = mysqlTable("test_results", {
  id: int("id").autoincrement().primaryKey(),
  analogId: int("analog_id").notNull(),
  testType: mysqlEnum("test_type", ["admet", "docking", "toxicity", "pkpd", "quantum"]).notNull(),
  testStatus: mysqlEnum("test_status", ["pending", "running", "completed", "failed"]).notNull(),
  results: text("results"), // JSON results
  errorMessage: text("error_message"),
  runBy: int("run_by").notNull(), // user ID
  createdAt: timestamp("created_at").defaultNow().notNull(),
  completedAt: timestamp("completed_at"),
});

export type TestResult = typeof testResults.$inferSelect;
export type InsertTestResult = typeof testResults.$inferInsert;

/**
 * Admin notifications table - stores alerts for high-confidence discoveries
 */
export const notifications = mysqlTable("notifications", {
  id: int("id").autoincrement().primaryKey(),
  userId: int("user_id").notNull(), // admin user to notify
  analogId: int("analog_id"),
  title: varchar("title", { length: 255 }).notNull(),
  message: text("message").notNull(),
  notificationType: mysqlEnum("notification_type", ["new-discovery", "high-confidence", "patent-alert", "system"]).notNull(),
  isRead: int("is_read").default(0).notNull(), // 0 = unread, 1 = read
  createdAt: timestamp("created_at").defaultNow().notNull(),
});

export type Notification = typeof notifications.$inferSelect;
export type InsertNotification = typeof notifications.$inferInsert;