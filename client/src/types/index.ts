export interface AnalogDiscovery {
  id: number;
  compoundId: string;
  compoundName: string;
  parentCompound: string;
  smiles: string;
  confidenceScore: number;
  similarityScore: number;
  safetyScore: number;
  efficacyScore: number;
  drugLikenessScore: number;
  patentStatus: "patent-free" | "patent-opportunity" | "patented" | "unknown";
  patentNumbers?: string;
  fdaStatus?: string;
  marketValue?: string;
  therapeuticPotential?: string;
  keyDifferences?: string;
  mechanismOfAction?: string;
  molecularWeight?: string;
  logP?: string;
  hBondDonors?: number;
  hBondAcceptors?: number;
  pubchemCid?: string;
  chemblId?: string;
  discoveredBy: string;
  discoveryMethod?: string;
  discoveredAt: Date;
  createdAt: Date;
  updatedAt: Date;
}

export interface TestResult {
  id: number;
  analogId: number;
  testType: "admet" | "docking" | "toxicity" | "pkpd" | "quantum";
  testStatus: "pending" | "running" | "completed" | "failed";
  results?: string;
  errorMessage?: string;
  runBy: number;
  createdAt: Date;
  completedAt?: Date;
}

export interface AnalyticsStats {
  totalDiscovered: number;
  highConfidenceCount: number;
  highConfidencePercentage: number;
  patentFreeCount: number;
  patentedCount: number;
  patentFreePercentage: number;
  averageConfidence: number;
}

export interface Notification {
  id: number;
  userId: number;
  analogId?: number;
  title: string;
  message: string;
  notificationType: "new-discovery" | "high-confidence" | "patent-alert" | "system";
  isRead: number;
  createdAt: Date;
}
