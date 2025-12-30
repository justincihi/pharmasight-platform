# Replit Merge Instructions - PharmaSight Dashboard

## Method 1: Fresh Clone (Recommended)

### Step 1: Create New Repl
1. Go to https://replit.com
2. Click **"+ Create Repl"**
3. Select **"Import from GitHub"**
4. Enter: `https://github.com/justincihi/pharmasight-platform`
5. Select branch: **`dashboard-manus`**
6. Name it: `pharmasight-dashboard`
7. Click **"Import from GitHub"**

### Step 2: Configure Environment
In Replit Secrets (🔒 icon in left sidebar), add:

```env
DATABASE_URL=mysql://user:password@host:port/database?ssl=true
JWT_SECRET=your-jwt-secret-here
GEMINI_API_KEY=your-gemini-key
ANTHROPIC_API_KEY=your-anthropic-key
SONAR_API_KEY=your-perplexity-key
PLATFORM_API_KEY=your-platform-api-key
```

### Step 3: Install & Run
```bash
pnpm install
pnpm db:push
pnpm dev
```

Your dashboard will be live at: `https://pharmasight-dashboard.your-username.repl.co`

---

## Method 2: Merge into Existing Repl

If you want to merge the dashboard into your existing `pharmasight-platform` Repl:

### Step 1: Pull Dashboard Branch
```bash
# In your Replit Shell
cd /home/runner/pharmasight-platform
git fetch origin dashboard-manus
git checkout dashboard-manus
```

### Step 2: Install Dashboard Dependencies
```bash
pnpm install
```

### Step 3: Set Up Environment Variables
Add the secrets listed above in Replit Secrets panel.

### Step 4: Push Database Schema
```bash
pnpm db:push
```

### Step 5: Run Dashboard
```bash
pnpm dev
```

---

## Method 3: Side-by-Side (Dashboard + Python Platform)

Keep both running in the same Repl:

### Directory Structure
```
pharmasight-platform/
├── src/                    # Python research engine
├── client/                 # Dashboard frontend
├── server/                 # Dashboard backend
├── drizzle/                # Database schema
└── package.json            # Dashboard dependencies
```

### Run Both Services

**Terminal 1 - Dashboard:**
```bash
pnpm dev
```

**Terminal 2 - Python Platform:**
```bash
python src/main.py
```

---

## Troubleshooting

### Port Conflicts
If port 3000 is taken, edit `server/_core/index.ts`:
```typescript
const PORT = process.env.PORT || 3001;
```

### Database Connection Issues
Ensure your DATABASE_URL includes SSL:
```
mysql://user:pass@host:port/db?ssl={"rejectUnauthorized":true}
```

### Missing Dependencies
```bash
rm -rf node_modules pnpm-lock.yaml
pnpm install
```

### Build Errors
```bash
pnpm build
```

---

## Video & PDF Upload Guide

### Adding Info Tab with Videos/PDFs

#### Step 1: Create Info Page Component

Create `client/src/pages/InfoHub.tsx`:

```typescript
import { Card } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";

export default function InfoHub() {
  return (
    <div className="container py-8">
      <h1 className="text-3xl font-bold mb-6">PharmaSight™ Info Hub</h1>
      
      <Tabs defaultValue="videos">
        <TabsList>
          <TabsTrigger value="videos">Videos</TabsTrigger>
          <TabsTrigger value="docs">Documentation</TabsTrigger>
          <TabsTrigger value="infographics">Infographics</TabsTrigger>
        </TabsList>

        <TabsContent value="videos">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* YouTube Embed */}
            <Card className="p-4">
              <h3 className="font-semibold mb-2">Platform Overview</h3>
              <iframe
                width="100%"
                height="315"
                src="https://www.youtube.com/embed/YOUR_VIDEO_ID"
                title="Platform Overview"
                frameBorder="0"
                allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
                allowFullScreen
              />
            </Card>

            {/* Vimeo Embed */}
            <Card className="p-4">
              <h3 className="font-semibold mb-2">Tutorial</h3>
              <iframe
                src="https://player.vimeo.com/video/YOUR_VIDEO_ID"
                width="100%"
                height="315"
                frameBorder="0"
                allow="autoplay; fullscreen; picture-in-picture"
                allowFullScreen
              />
            </Card>
          </div>
        </TabsContent>

        <TabsContent value="docs">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {/* PDF Links */}
            <Card className="p-4">
              <h3 className="font-semibold mb-2">User Guide</h3>
              <a 
                href="/docs/user-guide.pdf" 
                target="_blank"
                className="text-blue-600 hover:underline"
              >
                Download PDF →
              </a>
            </Card>
          </div>
        </TabsContent>

        <TabsContent value="infographics">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* Embedded PDF */}
            <Card className="p-4">
              <h3 className="font-semibold mb-2">Workflow Diagram</h3>
              <embed
                src="/infographics/workflow.pdf"
                type="application/pdf"
                width="100%"
                height="600px"
              />
            </Card>

            {/* Image Infographic */}
            <Card className="p-4">
              <h3 className="font-semibold mb-2">Process Overview</h3>
              <img 
                src="/infographics/process.png" 
                alt="Process Overview"
                className="w-full"
              />
            </Card>
          </div>
        </TabsContent>
      </Tabs>
    </div>
  );
}
```

#### Step 2: Add Route

In `client/src/App.tsx`, add:
```typescript
<Route path="/info" element={<InfoHub />} />
```

#### Step 3: Add Navigation Link

In `client/src/components/Navigation.tsx`, add to `navItems`:
```typescript
{ label: "Info Hub", path: "/info", icon: Info }
```

#### Step 4: Upload Files

**For PDFs and images:**
1. Place files in `client/public/docs/` or `client/public/infographics/`
2. Reference them with absolute paths: `/docs/filename.pdf`

**For videos:**
- **Option 1:** Upload to YouTube/Vimeo and embed
- **Option 2:** Use S3 storage (already configured in dashboard)

```typescript
// Upload video to S3
import { storagePut } from "@/server/storage";

const videoFile = await fetch(videoUrl).then(r => r.arrayBuffer());
const { url } = await storagePut(
  `videos/${filename}.mp4`,
  Buffer.from(videoFile),
  "video/mp4"
);
```

---

## File Upload API

To allow users to upload videos/PDFs through the dashboard:

### Backend (server/routers.ts)

```typescript
upload: protectedProcedure
  .input(z.object({
    filename: z.string(),
    contentType: z.string(),
    base64Data: z.string(),
  }))
  .mutation(async ({ input }) => {
    const { storagePut } = await import('./storage');
    
    const buffer = Buffer.from(input.base64Data, 'base64');
    const { url } = await storagePut(
      `uploads/${input.filename}`,
      buffer,
      input.contentType
    );
    
    return { url };
  }),
```

### Frontend Component

```typescript
function FileUpload() {
  const uploadMutation = trpc.upload.useMutation();

  const handleUpload = async (file: File) => {
    const reader = new FileReader();
    reader.onload = async (e) => {
      const base64 = e.target?.result?.toString().split(',')[1];
      const result = await uploadMutation.mutateAsync({
        filename: file.name,
        contentType: file.type,
        base64Data: base64!,
      });
      console.log('Uploaded:', result.url);
    };
    reader.readAsDataURL(file);
  };

  return (
    <input
      type="file"
      accept="video/*,application/pdf,image/*"
      onChange={(e) => e.target.files?.[0] && handleUpload(e.target.files[0])}
    />
  );
}
```

---

## Quick Reference

### Dashboard Structure
- **Frontend:** `client/src/` (React + Tailwind)
- **Backend:** `server/` (Express + tRPC)
- **Database:** `drizzle/schema.ts` (MySQL/TiDB)
- **Static Files:** `client/public/`

### Key URLs After Deployment
- Dashboard: `https://your-repl.repl.co`
- API: `https://your-repl.repl.co/api/trpc`
- Platform API: `https://your-repl.repl.co/api/platform`

### Common Commands
```bash
pnpm dev          # Start development server
pnpm build        # Build for production
pnpm db:push      # Update database schema
pnpm test         # Run tests
npx tsx trigger_research.ts  # Manual research run
```
