/**
 * Shared constants for docking and analysis
 */

export const COMMON_TARGETS = [
  { name: 'NMDA', description: 'NMDA Receptor (Glutamate)' },
  { name: '5HT2A', description: '5-HT2A Receptor (Serotonin)' },
  { name: 'D2', description: 'Dopamine D2 Receptor' },
  { name: 'D1', description: 'Dopamine D1 Receptor' },
  { name: 'GABA_A', description: 'GABA-A Receptor' },
  { name: 'mGluR5', description: 'mGluR5 Receptor' },
  { name: 'Sigma1', description: 'Sigma-1 Receptor' },
  { name: 'Custom', description: 'Custom Target' },
];

export const DEFAULT_DOCKING_PARAMS = {
  boxCenterX: 0,
  boxCenterY: 0,
  boxCenterZ: 0,
  boxSizeX: 25,
  boxSizeY: 25,
  boxSizeZ: 25,
  exhaustiveness: 8,
  numPoses: 9,
};

export const DOCKING_PRESETS = {
  NMDA: {
    name: 'NMDA Receptor',
    boxCenterX: 10.5,
    boxCenterY: 8.2,
    boxCenterZ: 12.1,
    boxSizeX: 20,
    boxSizeY: 20,
    boxSizeZ: 20,
    exhaustiveness: 8,
    numPoses: 9,
  },
  '5HT2A': {
    name: '5-HT2A Receptor',
    boxCenterX: 5.3,
    boxCenterY: 10.1,
    boxCenterZ: 8.9,
    boxSizeX: 22,
    boxSizeY: 22,
    boxSizeZ: 22,
    exhaustiveness: 8,
    numPoses: 9,
  },
  D2: {
    name: 'Dopamine D2',
    boxCenterX: 12.0,
    boxCenterY: 9.5,
    boxCenterZ: 11.3,
    boxSizeX: 25,
    boxSizeY: 25,
    boxSizeZ: 25,
    exhaustiveness: 8,
    numPoses: 9,
  },
};

export const SUPPORTED_FILE_TYPES = {
  PDB: { ext: '.pdb', mime: 'chemical/x-pdb' },
  SDF: { ext: '.sdf', mime: 'chemical/x-mdl-sdfile' },
  MOL: { ext: '.mol', mime: 'chemical/x-mdl-molfile' },
  VIDEO: { ext: ['.mp4', '.webm', '.mov'], mime: ['video/mp4', 'video/webm', 'video/quicktime'] },
  DIAGRAM: { ext: ['.png', '.jpg', '.svg', '.pdf'], mime: ['image/png', 'image/jpeg', 'image/svg+xml', 'application/pdf'] },
};

export const TOAST_MESSAGES = {
  UPLOAD_SUCCESS: (filename: string) => `✅ Successfully uploaded: ${filename}`,
  UPLOAD_FAILED: (error: string) => `❌ Upload failed: ${error}`,
  SAVE_SUCCESS: (item: string) => `✅ Successfully saved: ${item}`,
  SAVE_FAILED: (error: string) => `❌ Save failed: ${error}`,
  DELETE_SUCCESS: (item: string) => `✅ Successfully deleted: ${item}`,
  DELETE_FAILED: (error: string) => `❌ Delete failed: ${error}`,
  VALIDATION_ERROR: (message: string) => `⚠️ ${message}`,
  DOCKING_STARTED: (compound: string) => `🔬 Docking started for ${compound}`,
  DOCKING_COMPLETE: (affinity: number) => `✅ Docking complete! Affinity: ${affinity.toFixed(2)} kcal/mol`,
  DOCKING_FAILED: (error: string) => `❌ Docking failed: ${error}`,
  BATCH_SUBMITTED: (count: number) => `📦 Batch job submitted with ${count} compounds`,
  EXPORT_SUCCESS: (format: string, count: number) => `✅ Exported ${count} results to ${format.toUpperCase()}`,
  EXPORT_FAILED: (error: string) => `❌ Export failed: ${error}`,
};
