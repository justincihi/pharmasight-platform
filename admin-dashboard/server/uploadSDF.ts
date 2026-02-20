import { Router } from 'express';
import multer from 'multer';
import { join } from 'path';
import { randomBytes } from 'crypto';
import { tmpdir } from 'os';

const router = Router();

// Configure multer for file uploads
const storage = multer.diskStorage({
  destination: (req, file, cb) => {
    cb(null, join(tmpdir(), 'pharmasight-uploads'));
  },
  filename: (req, file, cb) => {
    const uniqueSuffix = randomBytes(8).toString('hex');
    cb(null, `${uniqueSuffix}-${file.originalname}`);
  },
});

const upload = multer({
  storage,
  limits: {
    fileSize: 10 * 1024 * 1024, // 10 MB
  },
  fileFilter: (req, file, cb) => {
    if (file.originalname.endsWith('.sdf')) {
      cb(null, true);
    } else {
      cb(new Error('Only SDF files are allowed'));
    }
  },
});

// Create upload directory if it doesn't exist
import { mkdirSync, existsSync } from 'fs';
const uploadDir = join(tmpdir(), 'pharmasight-uploads');
if (!existsSync(uploadDir)) {
  mkdirSync(uploadDir, { recursive: true });
}

router.post('/upload-sdf', upload.single('file'), (req, res) => {
  if (!req.file) {
    return res.status(400).json({ error: 'No file uploaded' });
  }

  res.json({
    success: true,
    filePath: req.file.path,
    filename: req.file.originalname,
    size: req.file.size,
  });
});

export default router;
