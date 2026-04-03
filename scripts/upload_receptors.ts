/**
 * Upload PDBQT receptor files to S3 and populate database
 * Run with: pnpm exec ts-node scripts/upload_receptors.ts
 */

import fs from 'fs';
import path from 'path';
import { getDb } from '../server/db';
import { receptorLibrary } from '../drizzle/schema';
import { storagePut } from '../server/storage';

interface ReceptorMetadata {
  name: string;
  pdb_id: string;
  description: string;
  pdb_file: string;
  pdbqt_file: string;
  source: string;
  downloaded_at: string;
}

async function main() {
  console.log('📤 Uploading Receptor Files to S3 and Database\n');
  console.log('='.repeat(60));

  try {
    // Read metadata
    const receptorsDir = path.join(process.cwd(), 'receptors');
    const metadataFile = path.join(receptorsDir, 'metadata.json');

    if (!fs.existsSync(metadataFile)) {
      console.error('❌ Metadata file not found. Run download_receptors.py first.');
      process.exit(1);
    }

    const metadata: ReceptorMetadata[] = JSON.parse(
      fs.readFileSync(metadataFile, 'utf-8')
    );

    console.log(`Found ${metadata.length} receptors to upload\n`);

    const db = await getDb();
    let uploadedCount = 0;

    for (const receptor of metadata) {
      try {
        console.log(`\n📥 Processing ${receptor.name}...`);

        const pdbqtPath = path.join(receptorsDir, receptor.pdbqt_file);

        if (!fs.existsSync(pdbqtPath)) {
          console.error(`  ❌ File not found: ${pdbqtPath}`);
          continue;
        }

        // Read file
        const fileContent = fs.readFileSync(pdbqtPath);
        const fileSize = fileContent.length;

        // Upload to S3
        const s3Key = `receptors/${receptor.name}/${receptor.pdbqt_file}`;
        console.log(`  Uploading to S3: ${s3Key}`);

        const { url: pdbqtUrl } = await storagePut(s3Key, fileContent, 'text/plain');

        console.log(`  ✅ Uploaded: ${pdbqtUrl}`);

        // Insert into database
        console.log(`  Inserting into database...`);

        await db
          .insert(receptorLibrary)
          .values({
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
          })
          .catch((err: any) => {
            if (err.message.includes('Duplicate')) {
              console.log(`  ℹ️  Already exists in database`);
            } else {
              throw err;
            }
          });

        console.log(`  ✅ Database record created`);
        uploadedCount++;
      } catch (error: any) {
        console.error(`  ❌ Error: ${error.message}`);
      }
    }

    console.log(`\n${'='.repeat(60)}`);
    console.log(`✅ Upload complete! ${uploadedCount}/${metadata.length} receptors uploaded\n`);

    process.exit(0);
  } catch (error: any) {
    console.error('❌ Fatal error:', error.message);
    process.exit(1);
  }
}

main();
