import { useEffect, useRef, useState } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Loader2, Download, RotateCcw, Maximize2 } from 'lucide-react';

// Declare 3Dmol global
declare const $3Dmol: any;

interface MolecularViewerProps {
  smiles?: string;
  sdfData?: string;
  pdbData?: string;
  title?: string;
  showControls?: boolean;
  height?: number;
  showInteractions?: boolean;
}

export default function MolecularViewer({
  smiles,
  sdfData,
  pdbData,
  title = '3D Molecular Structure',
  showControls = true,
  height = 500,
  showInteractions = false,
}: MolecularViewerProps) {
  const viewerContainerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [scriptLoaded, setScriptLoaded] = useState(false);

  // Load 3Dmol.js script
  useEffect(() => {
    if (typeof window !== 'undefined' && !(window as any).$3Dmol) {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
      script.async = true;
      script.onload = () => setScriptLoaded(true);
      script.onerror = () => setError('Failed to load 3Dmol.js library');
      document.head.appendChild(script);
    } else {
      setScriptLoaded(true);
    }
  }, []);

  // Initialize viewer
  useEffect(() => {
    if (!scriptLoaded || !viewerContainerRef.current) return;

    const initViewer = async () => {
      try {
        setIsLoading(true);
        setError(null);

        // Create viewer
        const viewer = (window as any).$3Dmol.createViewer(viewerContainerRef.current, {
          backgroundColor: 'white',
        });
        viewerRef.current = viewer;

        // Load structure
        if (pdbData) {
          viewer.addModel(pdbData, 'pdb');
        } else if (sdfData) {
          viewer.addModel(sdfData, 'sdf');
        } else if (smiles) {
          // For SMILES, we need to convert to 3D first
          // Use a simple approach: fetch from PubChem or generate via backend
          try {
            const response = await fetch(`/api/trpc/analog.convertSmilesToSDF?input=${encodeURIComponent(JSON.stringify({ smiles }))}`);
            if (response.ok) {
              const data = await response.json();
              if (data.result?.data?.sdf) {
                viewer.addModel(data.result.data.sdf, 'sdf');
              } else {
                throw new Error('No SDF data returned');
              }
            } else {
              throw new Error('Failed to convert SMILES');
            }
          } catch (err) {
            // Fallback: show a simple ball-and-stick representation
            console.warn('Could not convert SMILES to 3D, using 2D representation');
            setError('3D structure not available. Upload SDF file for 3D visualization.');
            setIsLoading(false);
            return;
          }
        }

        // Set style
        viewer.setStyle({}, {
          stick: { radius: 0.15 },
          sphere: { scale: 0.25 },
        });

        // Add surface if showing interactions
        if (showInteractions) {
          viewer.addSurface((window as any).$3Dmol.SurfaceType.VDW, {
            opacity: 0.7,
            color: 'white',
          });
        }

        viewer.zoomTo();
        viewer.render();
        setIsLoading(false);
      } catch (err: any) {
        console.error('Failed to initialize viewer:', err);
        setError(err.message || 'Failed to load molecular structure');
        setIsLoading(false);
      }
    };

    initViewer();

    return () => {
      if (viewerRef.current) {
        viewerRef.current.clear();
      }
    };
  }, [scriptLoaded, smiles, sdfData, pdbData, showInteractions]);

  const handleReset = () => {
    if (viewerRef.current) {
      viewerRef.current.zoomTo();
      viewerRef.current.render();
    }
  };

  const handleScreenshot = () => {
    if (viewerRef.current) {
      const imgData = viewerRef.current.pngURI();
      const link = document.createElement('a');
      link.href = imgData;
      link.download = 'molecular-structure.png';
      link.click();
    }
  };

  const handleFullscreen = () => {
    if (viewerContainerRef.current) {
      if (viewerContainerRef.current.requestFullscreen) {
        viewerContainerRef.current.requestFullscreen();
      }
    }
  };

  const toggleStyle = (style: 'stick' | 'sphere' | 'cartoon') => {
    if (!viewerRef.current) return;

    viewerRef.current.setStyle({}, {});

    switch (style) {
      case 'stick':
        viewerRef.current.setStyle({}, {
          stick: { radius: 0.15 },
          sphere: { scale: 0.25 },
        });
        break;
      case 'sphere':
        viewerRef.current.setStyle({}, {
          sphere: { scale: 0.4 },
        });
        break;
      case 'cartoon':
        viewerRef.current.setStyle({}, {
          cartoon: { color: 'spectrum' },
        });
        break;
    }

    viewerRef.current.render();
  };

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle>{title}</CardTitle>
          {showControls && !isLoading && !error && (
            <div className="flex gap-2">
              <Button
                variant="outline"
                size="sm"
                onClick={handleReset}
                title="Reset view"
              >
                <RotateCcw className="w-4 h-4" />
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={handleScreenshot}
                title="Download screenshot"
              >
                <Download className="w-4 h-4" />
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={handleFullscreen}
                title="Fullscreen"
              >
                <Maximize2 className="w-4 h-4" />
              </Button>
            </div>
          )}
        </div>
      </CardHeader>
      <CardContent>
        <div
          ref={viewerContainerRef}
          style={{
            width: '100%',
            height: `${height}px`,
            position: 'relative',
            border: '1px solid #e5e7eb',
            borderRadius: '8px',
            overflow: 'hidden',
          }}
        >
          {isLoading && (
            <div className="absolute inset-0 flex items-center justify-center bg-gray-50 z-10">
              <div className="text-center">
                <Loader2 className="w-8 h-8 animate-spin mx-auto text-blue-600 mb-2" />
                <p className="text-sm text-gray-600">Loading 3D structure...</p>
              </div>
            </div>
          )}

          {error && (
            <div className="absolute inset-0 flex items-center justify-center bg-yellow-50 z-10">
              <div className="text-center p-6">
                <p className="text-yellow-800 font-semibold mb-2">Structure Not Available</p>
                <p className="text-sm text-yellow-700">{error}</p>
              </div>
            </div>
          )}
        </div>

        {!isLoading && !error && (
          <div className="mt-4 space-y-3">
            <div className="flex gap-2">
              <Button
                variant="outline"
                size="sm"
                onClick={() => toggleStyle('stick')}
              >
                Stick
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={() => toggleStyle('sphere')}
              >
                Sphere
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={() => toggleStyle('cartoon')}
              >
                Cartoon
              </Button>
            </div>
            <div className="text-xs text-gray-500 space-y-1">
              <p>• Left click + drag: Rotate</p>
              <p>• Right click + drag: Pan</p>
              <p>• Scroll: Zoom</p>
            </div>
          </div>
        )}
      </CardContent>
    </Card>
  );
}
