#!/usr/bin/env node

/**
 * Test script to verify dock.py subprocess execution
 * Run with: node test-docking.mjs
 */

import { spawn } from 'child_process';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// Test data
const testSMILES = 'CC(=O)Nc1ccc(O)cc1'; // Acetaminophen (valid SMILES)
const testReceptorPath = path.join(__dirname, 'receptors', 'nmda_test.pdb');

console.log('🧪 Testing dock.py subprocess execution...\n');
console.log(`SMILES: ${testSMILES}`);
console.log(`Receptor: ${testReceptorPath}`);
console.log(`Python script: ${path.join(__dirname, 'scripts', 'dock.py')}\n`);

// Spawn Python process
const pythonProcess = spawn('python3', [
  path.join(__dirname, 'scripts', 'dock.py'),
  '--smiles', testSMILES,
  '--receptor', testReceptorPath,
  '--cx', '0',
  '--cy', '0',
  '--cz', '0',
  '--sx', '20',
  '--sy', '20',
  '--sz', '20',
  '--exhaustiveness', '8',
  '--num_poses', '5',
]);

let output = '';
let errorOutput = '';

pythonProcess.stdout.on('data', (data) => {
  output += data.toString();
  process.stdout.write(data);
});

pythonProcess.stderr.on('data', (data) => {
  errorOutput += data.toString();
  process.stderr.write(data);
});

pythonProcess.on('close', (code) => {
  console.log(`\n✅ Process exited with code ${code}\n`);

  if (code === 0) {
    try {
      // Extract JSON from output (may have Vina progress output before it)
      const jsonMatch = output.match(/\{"success".*\}/s);
      if (!jsonMatch) {
        throw new Error('No JSON found in output');
      }
      const result = JSON.parse(jsonMatch[0]);
      console.log('📊 Docking Result:');
      console.log(JSON.stringify(result, null, 2));

      if (result.success) {
        console.log('\n✅ DOCKING TEST PASSED!');
        console.log(`   Binding Affinity: ${result.binding_affinity} kcal/mol`);
        console.log(`   Number of Poses: ${result.num_poses}`);
      } else {
        console.log(`\n❌ DOCKING FAILED: ${result.error}`);
      }
    } catch (e) {
      console.log('❌ Failed to parse JSON output');
      console.log('Error:', e.message);
      console.log('Raw output:', output);
    }
  } else {
    console.log('❌ DOCKING TEST FAILED');
    console.log('Error output:', errorOutput);
  }
});

pythonProcess.on('error', (err) => {
  console.error('❌ Failed to spawn Python process:', err);
});
