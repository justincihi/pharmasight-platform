import { storagePut, storageGet } from "./storage";
import { getDb } from "./db";
import { pdbReceptors } from "../drizzle/schema";
import { eq, and } from "drizzle-orm";
import crypto from "crypto";
import fs from "fs/promises";
import path from "path";

// Receptors directory - relative to server root
const RECEPTORS_DIR = path.join(process.cwd(), "receptors");

export interface PDBFile {
  id: string;
  name: string;
  fileKey: string;
  url: string;
  uploadedBy: string;
  uploadedAt: Date;
  fileSize: number;
  targetName: string;
  description?: string;
}

/**
 * Upload a PDB file to S3 and store metadata in database
 */
export async function uploadPDBFile(
  pdbContent: Buffer | string,
  fileName: string,
  userId: string,
  targetName: string,
  description?: string
): Promise<PDBFile> {
  try {
    // Validate PDB format (basic check)
    const content = typeof pdbContent === "string" ? pdbContent : pdbContent.toString();
    if (!content.includes("ATOM") && !content.includes("HETATM")) {
      throw new Error("Invalid PDB file: must contain ATOM or HETATM records");
    }

    // Sanitize target name for filesystem use
    const sanitizedTargetName = targetName
      .toLowerCase()
      .replace(/\s+/g, "_")
      .replace(/[^a-z0-9_-]/g, "");

    // Ensure receptors directory exists
    await fs.mkdir(RECEPTORS_DIR, { recursive: true });

    // Write PDB file to disk
    const diskFilePath = path.join(RECEPTORS_DIR, `${sanitizedTargetName}.pdb`);
    await fs.writeFile(diskFilePath, content, "utf-8");
    console.log(`[PDB Storage] Wrote PDB file to: ${diskFilePath}`);

    // Generate unique file key for S3 backup
    const fileId = crypto.randomBytes(8).toString("hex");
    const fileKey = `pdb-receptors/${userId}/${fileId}-${fileName}`;
    const fileSize = Buffer.byteLength(content);

    // Upload to S3 as backup
    const { url } = await storagePut(fileKey, content, "chemical/x-pdb");

    // Store metadata in database
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    await db
      .insert(pdbReceptors)
      .values({
        id: fileId,
        name: fileName,
        fileKey: fileKey,
        url: url,
        uploadedBy: userId,
        uploadedAt: new Date(),
        fileSize: fileSize,
        targetName: targetName,
        description: description || null,
      });

    return {
      id: fileId,
      name: fileName,
      fileKey: fileKey,
      url: url,
      uploadedBy: userId,
      uploadedAt: new Date(),
      fileSize: fileSize,
      targetName: targetName,
      description: description,
    };
  } catch (error: any) {
    console.error("[PDB Upload Error]", error);
    throw new Error(`Failed to upload PDB file: ${error.message}`);
  }
}

/**
 * Get the path to a PDB receptor file on disk
 */
export function getPDBFilePath(targetName: string): string {
  const sanitizedTargetName = targetName
    .toLowerCase()
    .replace(/\s+/g, "_")
    .replace(/[^a-z0-9_-]/g, "");
  return path.join(RECEPTORS_DIR, `${sanitizedTargetName}.pdb`);
}

/**
 * Get all PDB files for a user
 */
export async function getUserPDBFiles(userId: string): Promise<PDBFile[]> {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const files = await db
      .select()
      .from(pdbReceptors)
      .where(eq(pdbReceptors.uploadedBy, userId))
      .orderBy(pdbReceptors.uploadedAt);

    return files.map((f) => ({
      id: f.id,
      name: f.name,
      fileKey: f.fileKey,
      url: f.url,
      uploadedBy: f.uploadedBy,
      uploadedAt: f.uploadedAt,
      fileSize: f.fileSize,
      targetName: f.targetName,
      description: f.description || undefined,
    }));
  } catch (error: any) {
    console.error("[PDB Fetch Error]", error);
    throw new Error(`Failed to fetch PDB files: ${error.message}`);
  }
}

/**
 * Get a specific PDB file by ID
 */
export async function getPDBFileById(fileId: string, userId: string): Promise<PDBFile | null> {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const file = await db
      .select()
      .from(pdbReceptors)
      .where(and(eq(pdbReceptors.id, fileId), eq(pdbReceptors.uploadedBy, userId)))
      .limit(1);

    if (!file || file.length === 0) return null;

    return {
      id: file[0].id,
      name: file[0].name,
      fileKey: file[0].fileKey,
      url: file[0].url,
      uploadedBy: file[0].uploadedBy,
      uploadedAt: file[0].uploadedAt,
      fileSize: file[0].fileSize,
      targetName: file[0].targetName,
      description: file[0].description || undefined,
    };
  } catch (error: any) {
    console.error("[PDB Fetch Error]", error);
    throw new Error(`Failed to fetch PDB file: ${error.message}`);
  }
}

/**
 * Get PDB file content from S3
 */
export async function getPDBFileContent(fileKey: string): Promise<string> {
  try {
    const { url } = await storageGet(fileKey);
    const response = await fetch(url);
    if (!response.ok) throw new Error(`Failed to fetch PDB file: ${response.statusText}`);
    return await response.text();
  } catch (error: any) {
    console.error("[PDB Content Fetch Error]", error);
    throw new Error(`Failed to fetch PDB content: ${error.message}`);
  }
}

/**
 * Delete a PDB file
 */
export async function deletePDBFile(fileId: string, userId: string): Promise<boolean> {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    // Get file to verify ownership and get fileKey
    const file = await db
      .select()
      .from(pdbReceptors)
      .where(and(eq(pdbReceptors.id, fileId), eq(pdbReceptors.uploadedBy, userId)))
      .limit(1);

    if (!file || file.length === 0) {
      throw new Error("PDB file not found or unauthorized");
    }

    // Delete from database
    await db.delete(pdbReceptors).where(eq(pdbReceptors.id, fileId));

    // Note: S3 deletion would require additional implementation
    // For now, we just remove the database record

    return true;
  } catch (error: any) {
    console.error("[PDB Delete Error]", error);
    throw new Error(`Failed to delete PDB file: ${error.message}`);
  }
}

/**
 * Get PDB files by target name
 */
export async function getPDBFilesByTarget(targetName: string): Promise<PDBFile[]> {
  try {
    const db = await getDb();
    if (!db) throw new Error("Database not available");

    const files = await db
      .select()
      .from(pdbReceptors)
      .where(eq(pdbReceptors.targetName, targetName))
      .orderBy(pdbReceptors.uploadedAt);

    return files.map((f) => ({
      id: f.id,
      name: f.name,
      fileKey: f.fileKey,
      url: f.url,
      uploadedBy: f.uploadedBy,
      uploadedAt: f.uploadedAt,
      fileSize: f.fileSize,
      targetName: f.targetName,
      description: f.description || undefined,
    }));
  } catch (error: any) {
    console.error("[PDB Fetch Error]", error);
    throw new Error(`Failed to fetch PDB files by target: ${error.message}`);
  }
}
