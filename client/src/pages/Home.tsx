import { useState } from "react";
import { useAuth } from "@/_core/hooks/useAuth";
import { Button } from "@/components/ui/button";
import { ChatbotInterface } from "@/components/ChatbotInterface";
import { getLoginUrl } from "@/const";
import { Link } from "wouter";
import {
  Beaker,
  Database,
  Zap,
  Microscope,
  Shield,
  TrendingUp,
  Loader2,
} from "lucide-react";

const INTEGRATIONS = [
  { name: "PubChem", icon: Database, description: "Compound database" },
  { name: "ChEMBL", icon: Microscope, description: "Bioactivity data" },
  { name: "RDKit", icon: Beaker, description: "Cheminformatics" },
  { name: "Biotransformer", icon: Zap, description: "Metabolite prediction" },
  { name: "AutoDock Vina", icon: Shield, description: "Molecular docking" },
  { name: "ADMET Predictor", icon: TrendingUp, description: "ADMET analysis" },
];

export default function Home() {
  const { user, loading, isAuthenticated } = useAuth();
  const [showChat, setShowChat] = useState(false);

  if (loading) {
    return (
      <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-blue-50 to-indigo-100">
        <Loader2 className="w-8 h-8 animate-spin text-blue-600" />
      </div>
    );
  }

  if (!isAuthenticated) {
    return (
      <div className="min-h-screen bg-gradient-to-br from-blue-50 via-white to-indigo-50">
        {/* Navigation */}
        <nav className="border-b border-gray-200 bg-white/80 backdrop-blur-sm sticky top-0 z-50">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4 flex items-center justify-between">
            <div className="flex items-center gap-3">
              <div className="w-10 h-10 bg-gradient-to-br from-blue-600 to-indigo-600 rounded-lg flex items-center justify-center">
                <Microscope className="w-6 h-6 text-white" />
              </div>
              <span className="text-xl font-bold text-gray-900">PharmaSight™</span>
            </div>
            <Button onClick={() => window.location.href = getLoginUrl()}>
              Sign In
            </Button>
          </div>
        </nav>

        {/* Hero Section */}
        <section className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-20">
          <div className="text-center mb-16">
            <h1 className="text-5xl md:text-6xl font-bold text-gray-900 mb-6">
              Autonomous Pharmaceutical
              <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-600 to-indigo-600">
                {" "}
                Analog Discovery
              </span>
            </h1>
            <p className="text-xl text-gray-600 mb-8 max-w-2xl mx-auto">
              Discover novel drug candidates with AI-powered cheminformatics.
              Analyze compounds, predict properties, and accelerate drug development.
            </p>
            <Button
              size="lg"
              onClick={() => window.location.href = getLoginUrl()}
              className="bg-gradient-to-r from-blue-600 to-indigo-600 hover:from-blue-700 hover:to-indigo-700"
            >
              Get Started
            </Button>
          </div>

          {/* Integration Badges */}
          <div className="bg-white rounded-xl border border-gray-200 p-8 mb-16">
            <h2 className="text-2xl font-bold text-gray-900 mb-8 text-center">
              Powered by Leading Cheminformatics Tools
            </h2>
            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-4">
              {INTEGRATIONS.map((integration) => {
                const Icon = integration.icon;
                return (
                  <div
                    key={integration.name}
                    className="flex flex-col items-center p-4 rounded-lg hover:bg-blue-50 transition-colors"
                  >
                    <Icon className="w-8 h-8 text-blue-600 mb-2" />
                    <p className="font-semibold text-sm text-gray-900 text-center">
                      {integration.name}
                    </p>
                    <p className="text-xs text-gray-500 text-center mt-1">
                      {integration.description}
                    </p>
                  </div>
                );
              })}
            </div>
          </div>

          {/* Features Grid */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8 mb-16">
            <div className="bg-white rounded-xl p-6 border border-gray-200 hover:shadow-lg transition-shadow">
              <div className="w-12 h-12 bg-blue-100 rounded-lg flex items-center justify-center mb-4">
                <Beaker className="w-6 h-6 text-blue-600" />
              </div>
              <h3 className="text-lg font-semibold text-gray-900 mb-2">
                Analog Discovery
              </h3>
              <p className="text-gray-600">
                AI-powered discovery of novel pharmaceutical analogs with confidence scoring
              </p>
            </div>

            <div className="bg-white rounded-xl p-6 border border-gray-200 hover:shadow-lg transition-shadow">
              <div className="w-12 h-12 bg-indigo-100 rounded-lg flex items-center justify-center mb-4">
                <Microscope className="w-6 h-6 text-indigo-600" />
              </div>
              <h3 className="text-lg font-semibold text-gray-900 mb-2">
                Cheminformatics Analysis
              </h3>
              <p className="text-gray-600">
                ADMET prediction, molecular docking, toxicity assessment, and more
              </p>
            </div>

            <div className="bg-white rounded-xl p-6 border border-gray-200 hover:shadow-lg transition-shadow">
              <div className="w-12 h-12 bg-purple-100 rounded-lg flex items-center justify-center mb-4">
                <TrendingUp className="w-6 h-6 text-purple-600" />
              </div>
              <h3 className="text-lg font-semibold text-gray-900 mb-2">
                Patent Intelligence
              </h3>
              <p className="text-gray-600">
                FDA Orange Book integration and patent status tracking for discovered compounds
              </p>
            </div>
          </div>

          {/* Video Placeholder */}
          <div className="bg-white rounded-xl border border-gray-200 p-8 mb-16">
            <h2 className="text-2xl font-bold text-gray-900 mb-6">
              See PharmaSight in Action
            </h2>
            <div className="aspect-video bg-gradient-to-br from-gray-100 to-gray-200 rounded-lg flex items-center justify-center">
              <div className="text-center">
                <div className="w-16 h-16 bg-blue-600 rounded-full flex items-center justify-center mx-auto mb-4">
                  <div className="w-0 h-0 border-l-8 border-l-transparent border-r-0 border-t-5 border-t-transparent border-b-5 border-b-transparent ml-1" />
                </div>
                <p className="text-gray-600 font-medium">
                  Video demo coming soon
                </p>
                <p className="text-gray-500 text-sm mt-2">
                  Sign in to access the full platform
                </p>
              </div>
            </div>
          </div>

          {/* CTA */}
          <div className="text-center">
            <p className="text-gray-600 mb-6">
              Ready to accelerate your drug discovery?
            </p>
            <Button
              size="lg"
              onClick={() => window.location.href = getLoginUrl()}
              className="bg-gradient-to-r from-blue-600 to-indigo-600 hover:from-blue-700 hover:to-indigo-700"
            >
              Sign In to PharmaSight
            </Button>
          </div>
        </section>
      </div>
    );
  }

  // Authenticated view
  return (
    <div className="min-h-screen bg-gray-50">
      {/* Navigation */}
      <nav className="border-b border-gray-200 bg-white sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 bg-gradient-to-br from-blue-600 to-indigo-600 rounded-lg flex items-center justify-center">
              <Microscope className="w-6 h-6 text-white" />
            </div>
            <span className="text-xl font-bold text-gray-900">PharmaSight™</span>
          </div>
          <div className="flex items-center gap-4">
            <span className="text-sm text-gray-600">Welcome, {user?.name}</span>
            <Link href="/admin/dashboard">
              <Button variant="outline" size="sm">
                Dashboard
              </Button>
            </Link>
          </div>
        </div>
      </nav>

      {/* Main Content */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Chatbot */}
          <div className="lg:col-span-2 h-[600px]">
            <ChatbotInterface />
          </div>

          {/* Quick Stats */}
          <div className="space-y-4">
            <div className="bg-white rounded-lg border border-gray-200 p-6">
              <h3 className="font-semibold text-gray-900 mb-4">Quick Access</h3>
              <div className="space-y-2">
                <Link href="/admin/dashboard">
                  <Button variant="outline" className="w-full justify-start">
                    View All Analogs
                  </Button>
                </Link>
                <Link href="/analytics">
                  <Button variant="outline" className="w-full justify-start">
                    Analytics Dashboard
                  </Button>
                </Link>
                <Link href="/testing">
                  <Button variant="outline" className="w-full justify-start">
                    Run Analysis
                  </Button>
                </Link>
              </div>
            </div>

            <div className="bg-gradient-to-br from-blue-50 to-indigo-50 rounded-lg border border-blue-200 p-6">
              <h3 className="font-semibold text-gray-900 mb-2">Latest Discoveries</h3>
              <p className="text-sm text-gray-600 mb-4">
                Check the analytics dashboard for the latest analog discoveries and their properties.
              </p>
              <Link href="/analytics">
                <Button size="sm" className="w-full">
                  View Analytics
                </Button>
              </Link>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
