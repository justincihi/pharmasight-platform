import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { getDb } from './db';
import { receptorLibrary } from '../drizzle/schema';
import { storagePut } from './storage';
import fs from 'fs';
import path from 'path';
import { eq } from 'drizzle-orm';

/**
 * Receptor Library Upload Test
 * Uploads PDBQT files to S3 and populates the database
 */

let db: any;

beforeAll(async () => {
  db = await getDb();
  expect(db).toBeDefined();
});

describe('Receptor Library - PDBQT Upload', () => {
  
  it('should upload PDBQT files to S3 and database', async () => {
    const receptorsDir = path.join(process.cwd(), 'receptors');
    const metadataFile = path.join(receptorsDir, 'metadata.json');

    // Check if metadata exists
    if (!fs.existsSync(metadataFile)) {
      console.log('⚠️  Metadata file not found, skipping upload test');
      return;
    }

    const metadata = JSON.parse(fs.readFileSync(metadataFile, 'utf-8'));
    console.log(`\n📤 Uploading ${metadata.length} receptors to S3 and database...\n`);

    let uploadedCount = 0;

    for (const receptor of metadata) {
      try {
        console.log(`  📥 Processing ${receptor.name}...`);

        const pdbqtPath = path.join(receptorsDir, receptor.pdbqt_file);

        if (!fs.existsSync(pdbqtPath)) {
          console.log(`    ⚠️  File not found: ${receptor.pdbqt_file}`);
          continue;
        }

        // Read file
        const fileContent = fs.readFileSync(pdbqtPath);
        const fileSize = fileContent.length;

        // Upload to S3
        const s3Key = `receptors/${receptor.name}/${receptor.pdbqt_file}`;
        const { url: pdbqtUrl } = await storagePut(s3Key, fileContent, 'text/plain');

        console.log(`    ✅ Uploaded to S3: ${pdbqtUrl}`);

        // Check if already exists
        const existing = await db
          .select()
          .from(receptorLibrary)
          .where(eq(receptorLibrary.targetName, receptor.name))
          .then((rows: any) => rows[0]);

        if (existing) {
          console.log(`    ℹ️  Already exists in database`);
          continue;
        }

        // Insert into database
        await db.insert(receptorLibrary).values({
          targetName: receptor.name,
          description: receptor.description,
          pdbId: receptor.pdb_id,
          pdbFileName: receptor.pdb_file,
          pdbqtFileName: receptor.pdbqt_file,
          pdbqtUrl: pdbqtUrl,
          fileSize: fileSize,
          category: 'psychiatric',
          source: receptor.source,
          uploadedBy: 'system',
          isActive: true,
        });

        console.log(`    ✅ Database record created`);
        uploadedCount++;
      } catch (error: any) {
        console.error(`    ❌ Error: ${error.message}`);
      }
    }

    console.log(`\n✅ Upload complete! ${uploadedCount}/${metadata.length} receptors uploaded\n`);

    // Verify receptors are in database
    const receptors = await db.select().from(receptorLibrary);
    expect(receptors.length).toBeGreaterThan(0);
  });

  it('should list all receptors', async () => {
    const receptors = await db.select().from(receptorLibrary);

    console.log(`\n📋 Receptor Library Summary:`);
    console.log(`   Total receptors: ${receptors.length}\n`);

    for (const receptor of receptors) {
      console.log(`   • ${receptor.targetName}`);
      console.log(`     PDB: ${receptor.pdbId}`);
      console.log(`     Description: ${receptor.description}`);
      console.log(`     File size: ${(receptor.fileSize / 1024).toFixed(2)} KB`);
    }

    expect(receptors.length).toBeGreaterThan(0);
  });

  it('should retrieve a specific receptor', async () => {
    const receptors = await db.select().from(receptorLibrary).limit(1);

    if (receptors.length === 0) {
      console.log('⚠️  No receptors in database');
      return;
    }

    const receptor = receptors[0];
    expect(receptor.targetName).toBeDefined();
    expect(receptor.pdbqtUrl).toBeDefined();
    expect(receptor.fileSize).toBeGreaterThan(0);

    console.log(`\n✅ Retrieved receptor: ${receptor.targetName}`);
    console.log(`   URL: ${receptor.pdbqtUrl}`);
    console.log(`   Size: ${(receptor.fileSize / 1024).toFixed(2)} KB`);
  });

  it('should get receptor statistics', async () => {
    const receptors = await db.select().from(receptorLibrary);

    const categories = new Set(receptors.map((r: any) => r.category));
    const totalSize = receptors.reduce((sum: number, r: any) => sum + (r.fileSize || 0), 0);

    const stats = {
      totalReceptors: receptors.length,
      categories: Array.from(categories),
      totalFileSize: totalSize,
      averageFileSize: receptors.length > 0 ? totalSize / receptors.length : 0,
    };

    console.log(`\n📊 Receptor Library Statistics:`);
    console.log(`   Total receptors: ${stats.totalReceptors}`);
    console.log(`   Categories: ${stats.categories.join(', ')}`);
    console.log(`   Total size: ${(stats.totalFileSize / 1024 / 1024).toFixed(2)} MB`);
    console.log(`   Average size: ${(stats.averageFileSize / 1024).toFixed(2)} KB`);

    expect(stats.totalReceptors).toBeGreaterThan(0);
  });
});
