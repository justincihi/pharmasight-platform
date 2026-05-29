import { int, mysqlEnum, mysqlTable, text, timestamp, varchar, tinyint, decimal, boolean, json } from "drizzle-orm/mysql-core";

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
  
  // Docking scores
  bindingAffinity: varchar("binding_affinity", { length: 64 }), // kcal/mol
  dockingScore: int("docking_score"), // 0-100 normalized score
  dockingTarget: varchar("docking_target", { length: 128 }), // e.g., "NMDA Receptor"
  
  // External database IDs
  pubchemCid: varchar("pubchem_cid", { length: 64 }),
  chemblId: varchar("chembl_id", { length: 64 }),
  
  // Discovery metadata
  discoveredBy: varchar("discovered_by", { length: 128 }).notNull(), // "autonomous-engine" or user ID
  discoveryMethod: varchar("discovery_method", { length: 128 }),
  discoveredAt: timestamp("discovered_at").defaultNow().notNull(),
  
  // Lead optimization tracking
  parentAnalogId: int("parent_analog_id"), // ID of the parent analog this was optimized from
  optimizationGeneration: int("optimization_generation").default(1), // Generation number (1 = original, 2 = first optimization, etc.)
  optimizationTarget: varchar("optimization_target", { length: 128 }), // What property was being optimized
  optimizationNotes: text("optimization_notes"), // Notes about the optimization
  
  // Advanced analysis results (JSON)
  toxicityProfile: text("toxicity_profile"), // JSON: hERG, hepatotoxicity, mutagenicity, carcinogenicity
  syntheticAccessibility: text("synthetic_accessibility"), // JSON: SA score, difficulty, estimated steps
  metabolites: text("metabolites"), // JSON: predicted metabolites
  
  // Approval workflow
  approvalStatus: mysqlEnum("approval_status", ["pending", "approved", "rejected"]).default("pending").notNull(),
  approvedBy: varchar("approved_by", { length: 128 }), // user openId who approved/rejected
  approvedAt: timestamp("approved_at"),
  
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

/**
 * Chat messages table - stores conversation history for context and recall
 */
export const chatMessages = mysqlTable("chat_messages", {
  id: int("id").autoincrement().primaryKey(),
  userId: int("user_id").notNull(), // user who sent the message
  role: mysqlEnum("role", ["user", "assistant", "system"]).notNull(),
  content: text("content").notNull(),
  provider: varchar("provider", { length: 64 }), // which LLM provider was used
  analogIds: text("analog_ids"), // JSON array of analog IDs mentioned
  createdAt: timestamp("created_at").defaultNow().notNull(),
});

export type ChatMessage = typeof chatMessages.$inferSelect;
export type InsertChatMessage = typeof chatMessages.$inferInsert;

/**
 * Docking queue table - manages automated docking jobs for multiple targets
 */
export const dockingQueue = mysqlTable("docking_queue", {
  id: int("id").autoincrement().primaryKey(),
  analogId: int("analog_id").notNull(),
  target: varchar("target", { length: 100 }).notNull(), // 'NMDA', '5-HT2A', 'D2'
  status: mysqlEnum("status", ["pending", "running", "completed", "failed"]).default("pending").notNull(),
  priority: int("priority").default(5).notNull(), // 1-10, higher = more important
  
  // Results
  bindingAffinity: varchar("binding_affinity", { length: 50 }),
  dockingScore: int("docking_score"),
  ligandPDB: text("ligand_pdb"),
  receptorPDB: text("receptor_pdb"),
  
  // Metadata
  errorMessage: text("error_message"),
  startedAt: timestamp("started_at"),
  completedAt: timestamp("completed_at"),
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type DockingQueueEntry = typeof dockingQueue.$inferSelect;
export type NewDockingQueueEntry = typeof dockingQueue.$inferInsert;

/**
 * Metabolites table - stores predicted metabolites for each analog
 */
export const metabolites = mysqlTable("metabolites", {
  id: int("id").autoincrement().primaryKey(),
  parentAnalogId: int("parent_analog_id").notNull(),
  smiles: text("smiles").notNull(),
  transformation: varchar("transformation", { length: 255 }).notNull(),
  phase: mysqlEnum("phase", ["Phase I", "Phase II"]).notNull(),
  enzyme: varchar("enzyme", { length: 64 }).notNull(),
  probability: varchar("probability", { length: 20 }).notNull(),
  molecularWeight: varchar("molecular_weight", { length: 20 }),
  logP: varchar("log_p", { length: 20 }),
  metabolicStabilityScore: int("metabolic_stability_score"),
  admetScore: int("admet_score"),
  dockingScore: int("docking_score"),
  createdAt: timestamp("created_at").defaultNow().notNull(),
});

export type Metabolite = typeof metabolites.$inferSelect;
export type InsertMetabolite = typeof metabolites.$inferInsert;

/**
 * Bookmarks table - allows users to save/bookmark important discoveries
 */
export const bookmarks = mysqlTable("bookmarks", {
  id: int("id").autoincrement().primaryKey(),
  userId: int("user_id").notNull(), // user who bookmarked
  analogId: int("analog_id"), // bookmarked analog discovery (optional)
  notificationId: int("notification_id"), // bookmarked notification (optional)
  title: varchar("title", { length: 255 }).notNull(),
  notes: text("notes"), // user's personal notes about this bookmark
  category: mysqlEnum("category", ["high-priority", "review-later", "promising", "archived"]).default("review-later").notNull(),
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type Bookmark = typeof bookmarks.$inferSelect;
export type InsertBookmark = typeof bookmarks.$inferInsert;


/**
 * PDB Receptor files table - stores uploaded protein receptor structures for docking
 */
export const pdbReceptors = mysqlTable("pdb_receptors", {
  id: varchar("id", { length: 64 }).primaryKey(), // Unique file ID
  name: varchar("name", { length: 255 }).notNull(), // Original filename
  fileKey: varchar("file_key", { length: 512 }).notNull(), // S3 storage key
  url: text("url").notNull(), // S3 presigned URL
  uploadedBy: varchar("uploaded_by", { length: 128 }).notNull(), // User ID who uploaded
  uploadedAt: timestamp("uploaded_at").defaultNow().notNull(),
  fileSize: int("file_size").notNull(), // File size in bytes
  targetName: varchar("target_name", { length: 255 }).notNull(), // e.g., "NMDA Receptor", "5HT2A"
  description: text("description"), // Optional description of the receptor
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type PDBReceptor = typeof pdbReceptors.$inferSelect;
export type InsertPDBReceptor = typeof pdbReceptors.$inferInsert;

/**
 * Docking parameters table - stores custom docking configurations for different targets
 */
export const dockingParameters = mysqlTable("docking_parameters", {
  id: int("id").autoincrement().primaryKey(),
  userId: varchar("user_id", { length: 128 }).notNull(), // User who created this configuration
  name: varchar("name", { length: 255 }).notNull(), // e.g., "NMDA Standard", "5HT2A High Exhaustiveness"
  targetName: varchar("target_name", { length: 255 }).notNull(), // Target protein name
  boxCenterX: varchar("box_center_x", { length: 64 }).notNull(), // X coordinate
  boxCenterY: varchar("box_center_y", { length: 64 }).notNull(), // Y coordinate
  boxCenterZ: varchar("box_center_z", { length: 64 }).notNull(), // Z coordinate
  boxSizeX: varchar("box_size_x", { length: 64 }).notNull(), // Box size X
  boxSizeY: varchar("box_size_y", { length: 64 }).notNull(), // Box size Y
  boxSizeZ: varchar("box_size_z", { length: 64 }).notNull(), // Box size Z
  exhaustiveness: int("exhaustiveness").notNull().default(8), // Vina exhaustiveness (1-32)
  numPoses: int("num_poses").notNull().default(9), // Number of poses to generate
  isDefault: tinyint("is_default").default(0), // Mark as default for target (0=false, 1=true)
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type DockingParameter = typeof dockingParameters.$inferSelect;
export type InsertDockingParameter = typeof dockingParameters.$inferInsert;

/**
 * Batch docking jobs table - tracks batch docking analysis jobs
 */
export const batchDockingJobs = mysqlTable("batch_docking_jobs", {
  id: varchar("id", { length: 64 }).primaryKey(), // Unique job ID
  userId: varchar("user_id", { length: 128 }).notNull(), // User who submitted the job
  jobName: varchar("job_name", { length: 255 }).notNull(), // User-provided job name
  status: mysqlEnum("status", ["pending", "running", "completed", "failed", "cancelled"]).default("pending").notNull(),
  totalCompounds: int("total_compounds").notNull(), // Total compounds to dock
  completedCompounds: int("completed_compounds").default(0).notNull(), // Compounds processed
  failedCompounds: int("failed_compounds").default(0).notNull(), // Compounds that failed
  targetName: varchar("target_name", { length: 255 }).notNull(), // Target protein
  parametersId: int("parameters_id"), // Reference to docking parameters used
  resultsSummary: text("results_summary"), // JSON summary of results
  createdAt: timestamp("created_at").defaultNow().notNull(),
  startedAt: timestamp("started_at"),
  completedAt: timestamp("completed_at"),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type BatchDockingJob = typeof batchDockingJobs.$inferSelect;
export type InsertBatchDockingJob = typeof batchDockingJobs.$inferInsert;

/**
 * Batch docking results table - stores individual results for each compound in a batch job
 */
export const batchDockingResults = mysqlTable("batch_docking_results", {
  id: int("id").autoincrement().primaryKey(),
  jobId: varchar("job_id", { length: 64 }).notNull(), // Reference to batch job
  analogId: int("analog_id").notNull(), // Reference to analog compound
  status: mysqlEnum("status", ["pending", "completed", "failed"]).default("pending").notNull(),
  bindingAffinity: varchar("binding_affinity", { length: 64 }), // kcal/mol
  dockingScore: int("docking_score"), // 0-100 normalized score
  numPoses: int("num_poses"), // Number of poses generated
  topPoses: text("top_poses"), // JSON array of top poses
  errorMessage: text("error_message"), // Error details if failed
  completedAt: timestamp("completed_at"),
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type BatchDockingResult = typeof batchDockingResults.$inferSelect;
export type InsertBatchDockingResult = typeof batchDockingResults.$inferInsert;

/**
 * Receptor library table - stores pre-built PDBQT receptor files for docking
 */
export const receptorLibrary = mysqlTable("receptor_library", {
  id: int("id").autoincrement().primaryKey(),
  targetName: varchar("target_name", { length: 128 }).notNull().unique(), // e.g., "NMDA", "5HT2A"
  description: text("description").notNull(), // Full description
  pdbId: varchar("pdb_id", { length: 16 }).notNull(), // PDB ID from RCSB
  pdbFileName: varchar("pdb_file_name", { length: 255 }).notNull(), // Original PDB file name
  pdbqtFileName: varchar("pdbqt_file_name", { length: 255 }).notNull(), // PDBQT file name
  pdbqtUrl: text("pdbqt_url").notNull(), // S3 URL to PDBQT file
  pdbUrl: text("pdb_url"), // S3 URL to original PDB file
  chain: varchar("chain", { length: 4 }).default("A"), // Primary chain to use
  ligandChain: varchar("ligand_chain", { length: 4 }), // Ligand chain if applicable
  fileSize: int("file_size"), // PDBQT file size in bytes
  resolution: varchar("resolution", { length: 16 }), // Crystal structure resolution
  experimentalMethod: varchar("experimental_method", { length: 64 }), // e.g., "X-RAY DIFFRACTION", "CRYO-EM"
  organism: varchar("organism", { length: 255 }), // Source organism
  
  // Docking parameters (default box settings for this target)
  defaultBoxCenterX: decimal("default_box_center_x", { precision: 10, scale: 3 }),
  defaultBoxCenterY: decimal("default_box_center_y", { precision: 10, scale: 3 }),
  defaultBoxCenterZ: decimal("default_box_center_z", { precision: 10, scale: 3 }),
  defaultBoxSizeX: decimal("default_box_size_x", { precision: 10, scale: 3 }),
  defaultBoxSizeY: decimal("default_box_size_y", { precision: 10, scale: 3 }),
  defaultBoxSizeZ: decimal("default_box_size_z", { precision: 10, scale: 3 }),
  defaultExhaustiveness: int("default_exhaustiveness").default(8),
  defaultNumPoses: int("default_num_poses").default(5),
  
  // Metadata
  category: varchar("category", { length: 64 }), // e.g., "psychiatric", "neurological", "cardiovascular"
  tags: text("tags"), // JSON array of tags for filtering
  notes: text("notes"), // Additional notes about the receptor
  source: varchar("source", { length: 128 }).default("RCSB PDB"), // Data source
  
  // Tracking
  uploadedBy: varchar("uploaded_by", { length: 128 }), // User ID who uploaded
  isActive: boolean("is_active").default(true),
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type ReceptorLibrary = typeof receptorLibrary.$inferSelect;
export type InsertReceptorLibrary = typeof receptorLibrary.$inferInsert;


/**
 * Cheminformatics results table - stores pipeline execution results
 */
export const cheminformaticsResults = mysqlTable("cheminformatics_results", {
  id: int("id").autoincrement().primaryKey(),
  userId: int("user_id").notNull().references(() => users.id),
  
  // Input parameters
  inputSmiles: text("input_smiles").notNull(),
  canonicalSmiles: text("canonical_smiles"),
  workflow: mysqlEnum("workflow", ["similarity", "brics", "validate", "full_pipeline"]).notNull(),
  threshold: decimal("threshold", { precision: 3, scale: 2 }).default("0.70"),
  maxHits: int("max_hits").default(25),
  
  // Results (stored as JSON for flexibility)
  results: json("results").$type<{
    success: boolean;
    error?: string;
    data?: {
      parent_smiles?: string;
      parent_cid?: number;
      parent_name?: string;
      hits?: Array<{
        cid: number;
        name: string;
        smiles: string;
        mw: number;
        tanimoto: number;
        patent_free?: boolean;
        patents?: string[];
      }>;
      total_hits?: number;
      analogs?: Array<{
        smiles: string;
        tanimoto: number;
        mw?: number;
      }>;
      total_generated?: number;
      master_list?: Array<{
        smiles: string;
        tanimoto: number;
        mw?: number;
        cid?: number;
        patent_free?: boolean;
        flag?: string;
      }>;
      total_candidates?: number;
      patent_free_count?: number;
      in_pubchem?: boolean;
      canonical_smiles?: string;
      mw?: number;
      num_atoms?: number;
      num_bonds?: number;
    };
  }>().notNull(),
  
  // Metadata
  executionTime: int("execution_time"), // milliseconds
  status: mysqlEnum("status", ["pending", "running", "completed", "failed"]).default("completed"),
  notes: text("notes"),
  
  // Tracking
  createdAt: timestamp("created_at").defaultNow().notNull(),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type CheminformaticsResult = typeof cheminformaticsResults.$inferSelect;
export type InsertCheminformaticsResult = typeof cheminformaticsResults.$inferInsert;


/**
 * Analysis Results table - persistent logging of all analysis results
 * Stores docking, toxicity, ADMET, and PK/PD analysis results with timestamps
 */
export const analysisResults = mysqlTable("analysis_results", {
  id: int("id").autoincrement().primaryKey(),
  analogId: int("analog_id"),
  analysisType: mysqlEnum("analysis_type", ["docking", "toxicity", "admet", "pkpd"]).notNull(),
  smiles: text("smiles").notNull(),
  target: varchar("target", { length: 128 }),
  result: json("result").$type<Record<string, unknown>>().notNull(),
  source: mysqlEnum("source", ["python", "api", "fallback"]).default("python").notNull(),
  executionTime: int("execution_time"), // milliseconds
  createdAt: timestamp("created_at").defaultNow().notNull(),
  createdBy: int("created_by"),
  updatedAt: timestamp("updated_at").defaultNow().onUpdateNow().notNull(),
});

export type AnalysisResult = typeof analysisResults.$inferSelect;
export type InsertAnalysisResult = typeof analysisResults.$inferInsert;
