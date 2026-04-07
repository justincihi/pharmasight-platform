import React, { useEffect, useRef } from 'react';
import { Card } from '@/components/ui/card';

interface DockingPose {
  mode: number;
  affinity: number;
  rmsd: number;
  pdbContent?: string;
}

interface DockingPoseViewer3DProps {
  pose: DockingPose;
  ligandPdb?: string;
  receptorPdb?: string;
  title?: string;
}

export const DockingPoseViewer3D: React.FC<DockingPoseViewer3DProps> = ({
  pose,
  ligandPdb,
  receptorPdb,
  title = `Docking Pose ${pose.mode}`,
}) => {
  const viewerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);

  useEffect(() => {
    if (!viewerRef.current) return;

    // Load 3Dmol.js library
    const script = document.createElement('script');
    script.src = 'https://3Dmol.org/build/3Dmol-min.js';
    script.async = true;

    script.onload = () => {
      // Initialize viewer
      const viewer = (window as any).$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'white',
        defaultcolors: (window as any).$3Dmol.rasmolColorMap,
      });

      viewerInstanceRef.current = viewer;

      // Load receptor if provided
      if (receptorPdb) {
        viewer.addModel(receptorPdb, 'pdb');
        viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
      }

      // Load ligand if provided
      if (ligandPdb) {
        viewer.addModel(ligandPdb, 'pdb');
        viewer.setStyle({ model: -1 }, { stick: { colorscheme: 'Jmol' } });
      }

      // Zoom to fit
      viewer.zoomTo();
      viewer.render();
    };

    document.head.appendChild(script);

    return () => {
      if (document.head.contains(script)) {
        document.head.removeChild(script);
      }
    };
  }, [receptorPdb, ligandPdb]);

  return (
    <Card className="p-4">
      <div className="mb-4">
        <h3 className="text-lg font-semibold">{title}</h3>
        <div className="grid grid-cols-2 gap-2 mt-2 text-sm">
          <div>
            <span className="text-gray-600">Binding Affinity:</span>
            <span className="ml-2 font-mono font-bold">{pose.affinity.toFixed(2)} kcal/mol</span>
          </div>
          <div>
            <span className="text-gray-600">RMSD:</span>
            <span className="ml-2 font-mono font-bold">{pose.rmsd.toFixed(2)} Å</span>
          </div>
        </div>
      </div>

      <div
        ref={viewerRef}
        style={{
          width: '100%',
          height: '500px',
          border: '1px solid #e5e7eb',
          borderRadius: '0.5rem',
          backgroundColor: '#f9fafb',
        }}
      />

      <div className="mt-4 text-xs text-gray-500">
        <p>💡 Tip: Use mouse to rotate, scroll to zoom, right-click to pan</p>
      </div>
    </Card>
  );
};

export default DockingPoseViewer3D;
