import { useEffect, useRef } from "react";
import { Card } from "@/components/ui/card";

interface MoleculeViewer3DProps {
  smiles: string;
  compoundName?: string;
  width?: number;
  height?: number;
}

export function MoleculeViewer3D({
  smiles,
  compoundName,
  width = 400,
  height = 300,
}: MoleculeViewer3DProps) {
  const viewerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);

  useEffect(() => {
    if (!viewerRef.current) return;

    // Dynamically import 3Dmol
    import("3dmol/build/3Dmol.js").then((module) => {
      const $3Dmol = (module as any).default || module;

      // Initialize viewer
      const viewer = $3Dmol.createViewer(viewerRef.current!, {
        backgroundColor: "white",
      });
      viewerInstanceRef.current = viewer;

      // Convert SMILES to 3D structure using RDKit-like approach
      // For now, we'll use a simple SDF format or fetch from PubChem
      fetchStructureFromPubChem(smiles)
        .then((sdf) => {
          viewer.addModel(sdf, "sdf");
          viewer.setStyle({}, { stick: { colorscheme: "Jmol" } });
          viewer.zoomTo();
          viewer.render();
        })
        .catch((error) => {
          console.error("Error loading 3D structure:", error);
          // Fallback: show a simple ball-and-stick model
          viewer.addModel(`${smiles}`, "smi");
          viewer.setStyle({}, { stick: { colorscheme: "Jmol" } });
          viewer.zoomTo();
          viewer.render();
        });
    });

    return () => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.clear();
      }
    };
  }, [smiles]);

  return (
    <Card className="p-4">
      {compoundName && (
        <h3 className="text-sm font-semibold text-gray-900 mb-2">
          {compoundName} - 3D Structure
        </h3>
      )}
      <div
        ref={viewerRef}
        style={{
          width: `${width}px`,
          height: `${height}px`,
          position: "relative",
        }}
        className="border border-gray-200 rounded-lg"
      />
      <p className="text-xs text-gray-500 mt-2">
        Drag to rotate • Scroll to zoom • Right-click to pan
      </p>
    </Card>
  );
}

// Helper function to fetch 3D structure from PubChem
async function fetchStructureFromPubChem(smiles: string): Promise<string> {
  try {
    // First, get the CID from SMILES
    const searchResponse = await fetch(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/cids/JSON`
    );
    const searchData = await searchResponse.json();
    const cid = searchData.IdentifierList.CID[0];

    // Then fetch the 3D SDF
    const sdfResponse = await fetch(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/SDF/?record_type=3d`
    );
    return await sdfResponse.text();
  } catch (error) {
    throw new Error("Failed to fetch 3D structure from PubChem");
  }
}
