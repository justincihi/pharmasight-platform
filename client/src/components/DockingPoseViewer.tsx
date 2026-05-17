import { useEffect, useRef, useState } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Loader2, Download, RotateCcw, Eye, EyeOff } from 'lucide-react';
import { Badge } from '@/components/ui/badge';
import DemoModeBadge, { DemoModeWarning } from './DemoModeBadge';

// Declare 3Dmol global
declare const $3Dmol: any;

interface DockingPoseViewerProps {
  ligandPDB?: string;
  receptorPDB?: string;
  bindingAffinity?: string;
  dockingScore?: number;
  target?: string;
  height?: number;
  source?: 'python' | 'fallback' | 'mock' | 'api';
  timestamp?: Date;
}

export default function DockingPoseViewer({
  ligandPDB,
  receptorPDB,
  bindingAffinity,
  dockingScore,
  target = 'NMDA Receptor',
  height = 600,
  source = 'fallback',
  timestamp,
}: DockingPoseViewerProps) {
  const viewerContainerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [scriptLoaded, setScriptLoaded] = useState(false);
  const [showProtein, setShowProtein] = useState(true);
  const [showLigand, setShowLigand] = useState(true);
  const [showInteractions, setShowInteractions] = useState(true);

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

        // Load receptor (protein)
        if (receptorPDB) {
          const proteinModel = viewer.addModel(receptorPDB, 'pdb');
          viewer.setStyle({ model: proteinModel }, {
            cartoon: { color: 'spectrum', opacity: 0.7 },
          });
        }

        // Load ligand
        if (ligandPDB) {
          const ligandModel = viewer.addModel(ligandPDB, 'pdb');
          viewer.setStyle({ model: ligandModel }, {
            stick: { radius: 0.2, colorscheme: 'greenCarbon' },
            sphere: { scale: 0.3, colorscheme: 'greenCarbon' },
          });

          // Add surface around ligand
          viewer.addSurface((window as any).$3Dmol.SurfaceType.VDW, {
            opacity: 0.5,
            color: 'lightblue',
          }, { model: ligandModel });
        }

        viewer.zoomTo();
        viewer.render();
        setIsLoading(false);
      } catch (err: any) {
        console.error('Failed to initialize viewer:', err);
        setError(err.message || 'Failed to load docking pose');
        setIsLoading(false);
      }
    };

    initViewer();

    return () => {
      if (viewerRef.current) {
        viewerRef.current.clear();
      }
    };
  }, [scriptLoaded, ligandPDB, receptorPDB]);

  // Toggle protein visibility
  useEffect(() => {
    if (!viewerRef.current) return;
    
    if (showProtein) {
      viewerRef.current.setStyle({ model: 0 }, {
        cartoon: { color: 'spectrum', opacity: 0.7 },
      });
    } else {
      viewerRef.current.setStyle({ model: 0 }, {});
    }
    viewerRef.current.render();
  }, [showProtein]);

  // Toggle ligand visibility
  useEffect(() => {
    if (!viewerRef.current) return;
    
    if (showLigand) {
      viewerRef.current.setStyle({ model: 1 }, {
        stick: { radius: 0.2, colorscheme: 'greenCarbon' },
        sphere: { scale: 0.3, colorscheme: 'greenCarbon' },
      });
    } else {
      viewerRef.current.setStyle({ model: 1 }, {});
    }
    viewerRef.current.render();
  }, [showLigand]);

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
      link.download = `docking-pose-${target}.png`;
      link.click();
    }
  };

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <div>
            <CardTitle>Docking Pose Visualization</CardTitle>
            <div className="flex gap-2 mt-2 flex-wrap">
              {bindingAffinity && (
                <Badge variant="outline">
                  Binding: {bindingAffinity} kcal/mol
                </Badge>
              )}
              {dockingScore && (
                <Badge variant={dockingScore >= 80 ? 'default' : 'secondary'}>
                  Score: {dockingScore}/100
                </Badge>
              )}
              <Badge variant="outline">{target}</Badge>
              <DemoModeBadge source={source} />
            </div>
          </div>
          {!isLoading && !error && (
            <div className="flex gap-2">
              <Button
                variant="outline"
                size="sm"
                onClick={() => setShowProtein(!showProtein)}
                title="Toggle protein"
              >
                {showProtein ? <Eye className="w-4 h-4" /> : <EyeOff className="w-4 h-4" />}
                <span className="ml-1 text-xs">Protein</span>
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={() => setShowLigand(!showLigand)}
                title="Toggle ligand"
              >
                {showLigand ? <Eye className="w-4 h-4" /> : <EyeOff className="w-4 h-4" />}
                <span className="ml-1 text-xs">Ligand</span>
              </Button>
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
                <p className="text-sm text-gray-600">Loading docking pose...</p>
              </div>
            </div>
          )}

          {error && (
            <div className="absolute inset-0 flex items-center justify-center bg-yellow-50 z-10">
              <div className="text-center p-6">
                <p className="text-yellow-800 font-semibold mb-2">Docking Data Not Available</p>
                <p className="text-sm text-yellow-700">{error}</p>
              </div>
            </div>
          )}
        </div>

        {!isLoading && !error && (
          <div className="mt-4 space-y-4">
            {source !== 'python' && source !== 'api' && (
              <DemoModeWarning source={source} />
            )}
            <div className="text-xs text-gray-500 space-y-1 mb-3">
              <p>• Left click + drag: Rotate</p>
              <p>• Right click + drag: Pan</p>
              <p>• Scroll: Zoom</p>
            </div>

            <div className="grid grid-cols-2 gap-4 text-sm">
              <div>
                <span className="font-semibold">Protein:</span>
                <span className="ml-2 text-muted-foreground">Cartoon (spectrum)</span>
              </div>
              <div>
                <span className="font-semibold">Ligand:</span>
                <span className="ml-2 text-muted-foreground">Stick + Sphere (green)</span>
              </div>
            </div>
            {timestamp && (
              <div className="text-xs text-gray-500 pt-2 border-t">
                <p>Generated: {timestamp.toLocaleString()}</p>
              </div>
            )}
          </div>
        )}
      </CardContent>
    </Card>
  );
}
