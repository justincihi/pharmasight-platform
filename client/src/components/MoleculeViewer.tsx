import { useEffect, useRef, useState } from "react";
import { GlassCard } from "./GlassCard";
import { Loader2 } from "lucide-react";

interface MoleculeViewerProps {
  smiles?: string;
  sdf?: string;
  pdb?: string;
  width?: string;
  height?: string;
  style?: "stick" | "sphere" | "line" | "cartoon";
  backgroundColor?: string;
}

export function MoleculeViewer({
  smiles,
  sdf,
  pdb,
  width = "100%",
  height = "400px",
  style = "stick",
  backgroundColor = "transparent",
}: MoleculeViewerProps) {
  const viewerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!viewerRef.current) return;

    // Dynamically import 3Dmol to avoid SSR issues
    import("3dmol").then(($3Dmol) => {
      try {
        const viewer = $3Dmol.createViewer(viewerRef.current!, {
          backgroundColor: backgroundColor,
        });

        if (smiles) {
          // Convert SMILES to 3D structure using RDKit or external service
          // For now, show a placeholder message
          setError("SMILES visualization requires backend conversion to 3D coordinates");
          setLoading(false);
          return;
        }

        if (sdf) {
          viewer.addModel(sdf, "sdf");
        } else if (pdb) {
          viewer.addModel(pdb, "pdb");
        } else {
          setError("No molecular data provided");
          setLoading(false);
          return;
        }

        // Apply visualization style
        const styleConfig: any = {};
        if (style === "stick") {
          viewer.setStyle({}, { stick: { radius: 0.15 } });
        } else if (style === "sphere") {
          viewer.setStyle({}, { sphere: { scale: 0.3 } });
        } else if (style === "line") {
          viewer.setStyle({}, { line: {} });
        } else if (style === "cartoon") {
          viewer.setStyle({}, { cartoon: { color: "spectrum" } });
        }

        viewer.zoomTo();
        viewer.render();
        setLoading(false);
      } catch (err: any) {
        setError(err.message || "Failed to render molecule");
        setLoading(false);
      }
    }).catch((err) => {
      setError("Failed to load 3Dmol.js library");
      setLoading(false);
    });
  }, [smiles, sdf, pdb, style, backgroundColor]);

  return (
    <GlassCard className="p-4">
      <div className="relative" style={{ width, height }}>
        {loading && (
          <div className="absolute inset-0 flex items-center justify-center">
            <Loader2 className="w-8 h-8 animate-spin text-blue-500" />
          </div>
        )}
        {error && (
          <div className="absolute inset-0 flex items-center justify-center">
            <p className="text-sm text-muted-foreground text-center px-4">
              {error}
            </p>
          </div>
        )}
        <div
          ref={viewerRef}
          className="w-full h-full rounded-lg"
          style={{ display: loading || error ? "none" : "block" }}
        />
      </div>
    </GlassCard>
  );
}
