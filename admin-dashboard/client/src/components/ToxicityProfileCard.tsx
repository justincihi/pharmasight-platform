import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { AlertTriangle, Heart, Activity, Dna, Shield } from "lucide-react";
import { cn } from "@/lib/utils";

interface ToxicityEndpoint {
  risk_score: number;
  risk_level: string;
  recommendation: string;
  prediction?: string;
}

interface ToxicityProfileData {
  hERG: ToxicityEndpoint & {
    logP: number;
    molecular_weight: number;
    tpsa: number;
    aromatic_rings: number;
    basic_groups: number;
    structural_alerts: number;
  };
  hepatotoxicity: ToxicityEndpoint & {
    logP: number;
    molecular_weight: number;
    structural_alerts: number;
  };
  mutagenicity: ToxicityEndpoint & {
    structural_alerts: number;
    alert_types: string[];
  };
  carcinogenicity: ToxicityEndpoint & {
    structural_alerts: number;
    alert_types: string[];
    aromatic_rings: number;
  };
}

interface ToxicityProfileCardProps {
  profile: ToxicityProfileData;
  compact?: boolean;
}

export function ToxicityProfileCard({ profile, compact = false }: ToxicityProfileCardProps) {
  const getRiskColor = (level: string) => {
    switch (level.toLowerCase()) {
      case 'low':
        return 'bg-green-100 text-green-800 border-green-200';
      case 'medium':
        return 'bg-yellow-100 text-yellow-800 border-yellow-200';
      case 'high':
        return 'bg-red-100 text-red-800 border-red-200';
      default:
        return 'bg-gray-100 text-gray-800 border-gray-200';
    }
  };

  const getRiskIcon = (endpoint: string) => {
    switch (endpoint) {
      case 'hERG':
        return <Heart className="h-4 w-4" />;
      case 'hepatotoxicity':
        return <Activity className="h-4 w-4" />;
      case 'mutagenicity':
        return <Dna className="h-4 w-4" />;
      case 'carcinogenicity':
        return <Shield className="h-4 w-4" />;
      default:
        return <AlertTriangle className="h-4 w-4" />;
    }
  };

  const endpoints = [
    { key: 'hERG', label: 'Cardiac (hERG)', data: profile.hERG, description: 'Risk of cardiac arrhythmias' },
    { key: 'hepatotoxicity', label: 'Liver Toxicity', data: profile.hepatotoxicity, description: 'Risk of liver damage' },
    { key: 'mutagenicity', label: 'Mutagenicity', data: profile.mutagenicity, description: 'DNA mutation potential' },
    { key: 'carcinogenicity', label: 'Carcinogenicity', data: profile.carcinogenicity, description: 'Cancer risk potential' },
  ];

  if (compact) {
    return (
      <Card>
        <CardHeader className="pb-3">
          <CardTitle className="text-base flex items-center gap-2">
            <AlertTriangle className="h-4 w-4" />
            Toxicity Profile
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 gap-2">
            {endpoints.map((endpoint) => (
              <div key={endpoint.key} className="flex items-center justify-between p-2 rounded-lg border">
                <div className="flex items-center gap-2">
                  {getRiskIcon(endpoint.key)}
                  <span className="text-sm font-medium">{endpoint.label}</span>
                </div>
                <Badge className={cn("text-xs", getRiskColor(endpoint.data.risk_level))}>
                  {endpoint.data.risk_level}
                </Badge>
              </div>
            ))}
          </div>
        </CardContent>
      </Card>
    );
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <AlertTriangle className="h-5 w-5" />
          Comprehensive Toxicity Profile
        </CardTitle>
        <CardDescription>
          Predicted toxicity endpoints based on structural analysis
        </CardDescription>
      </CardHeader>
      <CardContent>
        <div className="space-y-4">
          {endpoints.map((endpoint) => (
            <div key={endpoint.key} className="border rounded-lg p-4">
              <div className="flex items-start justify-between mb-2">
                <div className="flex items-center gap-2">
                  {getRiskIcon(endpoint.key)}
                  <div>
                    <h4 className="font-semibold">{endpoint.label}</h4>
                    <p className="text-xs text-muted-foreground">{endpoint.description}</p>
                  </div>
                </div>
                <Badge className={cn("ml-2", getRiskColor(endpoint.data.risk_level))}>
                  {endpoint.data.risk_level}
                </Badge>
              </div>
              
              <div className="mt-3 space-y-2">
                <div className="flex items-center gap-2">
                  <div className="flex-1 bg-gray-200 rounded-full h-2">
                    <div
                      className={cn(
                        "h-2 rounded-full transition-all",
                        endpoint.data.risk_level === 'Low' && "bg-green-500",
                        endpoint.data.risk_level === 'Medium' && "bg-yellow-500",
                        endpoint.data.risk_level === 'High' && "bg-red-500"
                      )}
                      style={{ width: `${endpoint.data.risk_score}%` }}
                    />
                  </div>
                  <span className="text-sm font-medium w-12 text-right">
                    {endpoint.data.risk_score}
                  </span>
                </div>
                
                <p className="text-sm text-muted-foreground">
                  {endpoint.data.recommendation}
                </p>
                
                {endpoint.data.prediction && (
                  <p className="text-sm">
                    <span className="font-medium">Prediction:</span> {endpoint.data.prediction}
                  </p>
                )}
                
                {('structural_alerts' in endpoint.data) && endpoint.data.structural_alerts > 0 && (
                  <div className="flex items-center gap-2 text-sm">
                    <AlertTriangle className="h-3 w-3 text-yellow-600" />
                    <span className="text-yellow-600">
                      {endpoint.data.structural_alerts} structural alert{endpoint.data.structural_alerts > 1 ? 's' : ''} detected
                    </span>
                  </div>
                )}
              </div>
            </div>
          ))}
        </div>
      </CardContent>
    </Card>
  );
}

export default ToxicityProfileCard;
