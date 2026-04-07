#!/usr/bin/env node

/**
 * Upload PDBQT receptor files to S3 and populate database
 * Run with: node scripts/upload_receptors.mjs
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { spawn } from 'child_process';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const projectRoot = path.dirname(__dirname);
const receptorsDir = path.join(projectRoot, 'receptors');
const metadataFile = path.join(receptorsDir, 'metadata.json');

console.log('📤 Uploading Receptor Files to S3\n');
console.log('=' .repeat(60));

// Read metadata
if (!fs.existsSync(metadataFile)) {
  console.error('❌ Metadata file not found. Run download_receptors.py first.');
  process.exit(1);
}

const metadata = JSON.parse(fs.readFileSync(metadataFile, 'utf-8'));

console.log(`Found ${metadata.length} receptors to upload\n`);

// Create a Node.js script that will use the server's storage helpers
const uploadScript = `
import { getDb } from './server/db.js';
import { receptorLibrary } from './drizzle/schema.js';
import { storagePut } from './server/storage.js';
import fs from 'fs';
import path from 'path';

const metadata = ${JSON.stringify(metadata)};
const receptorsDir = '${receptorsDir}';

async function main() {
  const db = await getDb();
  
  for (const receptor of metadata) {
    try {
      console.log(\`\\n📥 Uploading \${receptor.name}...\`);
      
      const pdbqtPath = path.join(receptorsDir, receptor.pdbqt_file);
      
      if (!fs.existsSync(pdbqtPath)) {
        console.error(\`  ❌ File not found: \${pdbqtPath}\`);
        continue;
      }
      
      // Read file
      const fileContent = fs.readFileSync(pdbqtPath);
      const fileSize = fileContent.length;
      
      // Upload to S3
      const s3Key = \`receptors/\${receptor.name}/\${receptor.pdbqt_file}\`;
      console.log(\`  Uploading to S3: \${s3Key}\`);
      
      const { url: pdbqtUrl } = await storagePut(s3Key, fileContent, 'text/plain');
      
      console.log(\`  ✅ Uploaded: \${pdbqtUrl}\`);
      
      // Insert into database
      console.log(\`  Inserting into database...\`);
      
      await db.insert(receptorLibrary).values({
        targetName: receptor.name,
        description: receptor.description,
        pdbId: receptor.pdb_id,
        pdbFileName: receptor.pdb_file,
        pdbqtFileName: receptor.pdbqt_file,
        pdbqtUrl: pdbqtUrl,
        fileSize: fileSize,
        category: 'psychiatric',
        source: 'RCSB PDB',
        uploadedBy: 'system',
        isActive: true,
      }).catch(err => {
        if (err.message.includes('Duplicate')) {
          console.log(\`  ℹ️  Already exists in database\`);
        } else {
          throw err;
        }
      });
      
      console.log(\`  ✅ Database record created\`);
      
    } catch (error) {
      console.error(\`  ❌ Error: \${error.message}\`);
    }
  }
  
  console.log(\`\\n✅ Upload complete!\`);
  process.exit(0);
}

main().catch(err => {
  console.error('❌ Fatal error:', err);
  process.exit(1);
});
`;

// Write and execute the upload script
const uploadScriptPath = path.join(projectRoot, 'scripts', 'upload_receptors_internal.mjs');
fs.writeFileSync(uploadScriptPath, uploadScript);

console.log('Executing upload script...\n');

// Run the script using pnpm exec
const proc = spawn('pnpm', ['exec', 'node', uploadScriptPath], {
  cwd: projectRoot,
  stdio: 'inherit',
});

proc.on('exit', (code) => {
  // Clean up
  fs.unlinkSync(uploadScriptPath);
  process.exit(code);
});
