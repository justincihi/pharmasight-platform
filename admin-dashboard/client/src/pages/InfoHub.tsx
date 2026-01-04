import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { trpc } from "@/lib/trpc";
import { FileVideo, FileText, Upload, Trash2, ExternalLink, Play } from "lucide-react";
import { useState } from "react";
import { toast } from "sonner";

export default function InfoHub() {
  const [uploadingVideo, setUploadingVideo] = useState(false);
  const [uploadingPdf, setUploadingPdf] = useState(false);
  
  const { data: videos, refetch: refetchVideos } = trpc.infohub.listVideos.useQuery();
  const { data: pdfs, refetch: refetchPdfs } = trpc.infohub.listPdfs.useQuery();
  
  const uploadVideoMutation = trpc.infohub.uploadVideo.useMutation({
    onSuccess: () => {
      toast.success("Video uploaded successfully");
      refetchVideos();
      setUploadingVideo(false);
    },
    onError: (error) => {
      toast.error(`Upload failed: ${error.message}`);
      setUploadingVideo(false);
    },
  });

  const uploadPdfMutation = trpc.infohub.uploadPdf.useMutation({
    onSuccess: () => {
      toast.success("PDF uploaded successfully");
      refetchPdfs();
      setUploadingPdf(false);
    },
    onError: (error) => {
      toast.error(`Upload failed: ${error.message}`);
      setUploadingPdf(false);
    },
  });

  const deleteVideoMutation = trpc.infohub.deleteVideo.useMutation({
    onSuccess: () => {
      toast.success("Video deleted");
      refetchVideos();
    },
  });

  const deletePdfMutation = trpc.infohub.deletePdf.useMutation({
    onSuccess: () => {
      toast.success("PDF deleted");
      refetchPdfs();
    },
  });

  const handleVideoUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    if (!file.type.startsWith("video/")) {
      toast.error("Please upload a video file");
      return;
    }

    setUploadingVideo(true);
    
    // Convert file to base64
    const reader = new FileReader();
    reader.onload = () => {
      const base64 = reader.result as string;
      uploadVideoMutation.mutate({
        title: file.name.replace(/\.[^/.]+$/, ""),
        url: base64,
        description: "",
      });
    };
    reader.readAsDataURL(file);
  };

  const handlePdfUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    if (file.type !== "application/pdf") {
      toast.error("Please upload a PDF file");
      return;
    }

    setUploadingPdf(true);
    
    // Convert file to base64
    const reader = new FileReader();
    reader.onload = () => {
      const base64 = reader.result as string;
      uploadPdfMutation.mutate({
        title: file.name.replace(/\.[^/.]+$/, ""),
        url: base64,
        description: "",
      });
    };
    reader.readAsDataURL(file);
  };

  return (
    <div className="container mx-auto py-8 space-y-6">
      {/* Header */}
      <div>
        <h1 className="text-3xl font-bold">Info Hub</h1>
        <p className="text-muted-foreground mt-1">
          Platform demos, documentation, and educational resources
        </p>
      </div>

      <Tabs defaultValue="videos" className="w-full">
        <TabsList className="grid w-full grid-cols-2">
          <TabsTrigger value="videos">
            <FileVideo className="h-4 w-4 mr-2" />
            Videos
          </TabsTrigger>
          <TabsTrigger value="pdfs">
            <FileText className="h-4 w-4 mr-2" />
            PDFs
          </TabsTrigger>
        </TabsList>

        {/* Videos Tab */}
        <TabsContent value="videos" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Upload Video</CardTitle>
              <CardDescription>
                Add platform demos, tutorials, or educational videos
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="space-y-4">
                <div className="grid w-full max-w-sm items-center gap-1.5">
                  <Label htmlFor="video-upload">Video File</Label>
                  <Input
                    id="video-upload"
                    type="file"
                    accept="video/*"
                    onChange={handleVideoUpload}
                    disabled={uploadingVideo}
                  />
                </div>
                {uploadingVideo && (
                  <p className="text-sm text-muted-foreground">Uploading...</p>
                )}
              </div>
            </CardContent>
          </Card>

          {/* Video List */}
          <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
            {videos?.map((video: any) => (
              <Card key={video.id}>
                <CardHeader>
                  <div className="flex items-start justify-between">
                    <div className="space-y-1">
                      <CardTitle className="text-base">{video.title}</CardTitle>
                      {video.description && (
                        <CardDescription className="text-sm">
                          {video.description}
                        </CardDescription>
                      )}
                    </div>
                    <Button
                      variant="ghost"
                      size="icon"
                      onClick={() => deleteVideoMutation.mutate({ id: video.id })}
                    >
                      <Trash2 className="h-4 w-4" />
                    </Button>
                  </div>
                </CardHeader>
                <CardContent>
                  <div className="space-y-2">
                    <video
                      src={video.url}
                      controls
                      className="w-full rounded-md"
                      style={{ maxHeight: "200px" }}
                    />
                    <div className="flex gap-2">
                      <Button
                        variant="outline"
                        size="sm"
                        className="w-full"
                        onClick={() => window.open(video.url, "_blank")}
                      >
                        <ExternalLink className="h-3 w-3 mr-1" />
                        Open
                      </Button>
                    </div>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>

          {videos?.length === 0 && (
            <Card>
              <CardContent className="flex flex-col items-center justify-center py-12">
                <FileVideo className="h-12 w-12 text-muted-foreground mb-4" />
                <p className="text-muted-foreground">No videos uploaded yet</p>
              </CardContent>
            </Card>
          )}
        </TabsContent>

        {/* PDFs Tab */}
        <TabsContent value="pdfs" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Upload PDF</CardTitle>
              <CardDescription>
                Add documentation, research papers, or infographics
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="space-y-4">
                <div className="grid w-full max-w-sm items-center gap-1.5">
                  <Label htmlFor="pdf-upload">PDF File</Label>
                  <Input
                    id="pdf-upload"
                    type="file"
                    accept="application/pdf"
                    onChange={handlePdfUpload}
                    disabled={uploadingPdf}
                  />
                </div>
                {uploadingPdf && (
                  <p className="text-sm text-muted-foreground">Uploading...</p>
                )}
              </div>
            </CardContent>
          </Card>

          {/* PDF List */}
          <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
            {pdfs?.map((pdf: any) => (
              <Card key={pdf.id}>
                <CardHeader>
                  <div className="flex items-start justify-between">
                    <div className="space-y-1">
                      <CardTitle className="text-base">{pdf.title}</CardTitle>
                      {pdf.description && (
                        <CardDescription className="text-sm">
                          {pdf.description}
                        </CardDescription>
                      )}
                    </div>
                    <Button
                      variant="ghost"
                      size="icon"
                      onClick={() => deletePdfMutation.mutate({ id: pdf.id })}
                    >
                      <Trash2 className="h-4 w-4" />
                    </Button>
                  </div>
                </CardHeader>
                <CardContent>
                  <div className="space-y-2">
                    <div className="bg-muted rounded-md p-8 flex items-center justify-center">
                      <FileText className="h-16 w-16 text-muted-foreground" />
                    </div>
                    <Button
                      variant="outline"
                      size="sm"
                      className="w-full"
                      onClick={() => window.open(pdf.url, "_blank")}
                    >
                      <ExternalLink className="h-3 w-3 mr-1" />
                      View PDF
                    </Button>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>

          {pdfs?.length === 0 && (
            <Card>
              <CardContent className="flex flex-col items-center justify-center py-12">
                <FileText className="h-12 w-12 text-muted-foreground mb-4" />
                <p className="text-muted-foreground">No PDFs uploaded yet</p>
              </CardContent>
            </Card>
          )}
        </TabsContent>
      </Tabs>
    </div>
  );
}
