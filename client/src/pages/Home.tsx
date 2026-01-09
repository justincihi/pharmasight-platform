import { useState } from "react";
import { motion } from "framer-motion";
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
  ArrowRight,
  Sparkles,
} from "lucide-react";

const INTEGRATIONS = [
  { name: "PubChem", icon: Database, description: "Compound database", color: "from-blue-500 to-cyan-500" },
  { name: "ChEMBL", icon: Microscope, description: "Bioactivity data", color: "from-purple-500 to-pink-500" },
  { name: "RDKit", icon: Beaker, description: "Cheminformatics", color: "from-green-500 to-emerald-500" },
  { name: "Biotransformer", icon: Zap, description: "Metabolite prediction", color: "from-yellow-500 to-orange-500" },
  { name: "AutoDock Vina", icon: Shield, description: "Molecular docking", color: "from-red-500 to-rose-500" },
  { name: "ADMET Predictor", icon: TrendingUp, description: "ADMET analysis", color: "from-indigo-500 to-blue-500" },
];

const FEATURES = [
  {
    icon: Beaker,
    title: "Analog Discovery",
    description: "AI-powered discovery of novel pharmaceutical analogs with confidence scoring and structural analysis",
    gradient: "from-blue-500 to-cyan-500",
  },
  {
    icon: Microscope,
    title: "Cheminformatics Analysis",
    description: "ADMET prediction, molecular docking, toxicity assessment, and comprehensive property profiling",
    gradient: "from-purple-500 to-pink-500",
  },
  {
    icon: TrendingUp,
    title: "Patent Intelligence",
    description: "FDA Orange Book integration and patent status tracking for discovered compounds with real-time updates",
    gradient: "from-indigo-500 to-blue-500",
  },
];

export default function Home() {
  const { user, loading, isAuthenticated } = useAuth();

  if (loading) {
    return (
      <div className="min-h-screen flex items-center justify-center bg-gradient-animated">
        <Loader2 className="w-8 h-8 animate-spin text-white" />
      </div>
    );
  }

  if (!isAuthenticated) {
    return (
      <div className="min-h-screen bg-background relative overflow-hidden">
        {/* Animated Background */}
        <div className="absolute inset-0 bg-gradient-animated opacity-20" />
        <div className="absolute inset-0 bg-[radial-gradient(circle_at_30%_20%,oklch(0.7_0.25_250/20%),transparent_50%)]" />
        <div className="absolute inset-0 bg-[radial-gradient(circle_at_70%_80%,oklch(0.65_0.22_280/20%),transparent_50%)]" />

        {/* Navigation */}
        <motion.nav
          initial={{ y: -100 }}
          animate={{ y: 0 }}
          className="glass sticky top-0 z-50 border-b border-glass-border"
        >
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4 flex items-center justify-between">
            <motion.div
              className="flex items-center gap-3"
              whileHover={{ scale: 1.05 }}
            >
              <div className="w-10 h-10 bg-gradient-to-br from-blue-600 to-indigo-600 rounded-xl flex items-center justify-center shadow-lg">
                <Microscope className="w-6 h-6 text-white" />
              </div>
              <span className="text-xl font-bold text-foreground">PharmaSight™</span>
            </motion.div>
            <Button
              onClick={() => window.location.href = getLoginUrl()}
              className="bg-gradient-to-r from-blue-600 to-indigo-600 hover:from-blue-700 hover:to-indigo-700"
            >
              Sign In
            </Button>
          </div>
        </motion.nav>

        {/* Hero Section */}
        <section className="relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-20 md:py-32">
          <div className="text-center mb-20">
            <motion.div
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8 }}
              className="inline-flex items-center gap-2 px-4 py-2 rounded-full glass border border-glass-border mb-8"
            >
              <Sparkles className="w-4 h-4 text-blue-500" />
              <span className="text-sm font-medium text-muted-foreground">
                AI-Powered Drug Discovery Platform
              </span>
            </motion.div>

            <motion.h1
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.1 }}
              className="text-5xl md:text-7xl font-bold text-foreground mb-6 leading-tight"
            >
              Autonomous Pharmaceutical
              <br />
              <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-600 via-indigo-600 to-purple-600">
                Analog Discovery
              </span>
            </motion.h1>

            <motion.p
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.2 }}
              className="text-xl text-muted-foreground mb-10 max-w-3xl mx-auto leading-relaxed"
            >
              Discover novel drug candidates with AI-powered cheminformatics.
              Analyze compounds, predict properties, and accelerate drug development
              with cutting-edge computational tools.
            </motion.p>

            <motion.div
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.8, delay: 0.3 }}
              className="flex flex-col sm:flex-row items-center justify-center gap-4"
            >
              <Button
                size="lg"
                onClick={() => window.location.href = getLoginUrl()}
                className="bg-gradient-to-r from-blue-600 to-indigo-600 hover:from-blue-700 hover:to-indigo-700 text-lg px-8 py-6 group"
              >
                Get Started
                <ArrowRight className="w-5 h-5 ml-2 group-hover:translate-x-1 transition-transform" />
              </Button>
              <Button
                size="lg"
                variant="outline"
                className="glass border-glass-border text-lg px-8 py-6"
              >
                View Demo
              </Button>
            </motion.div>
          </div>

          {/* Integration Badges */}
          <motion.div
            initial={{ opacity: 0, y: 40 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.4 }}
            className="glass-strong rounded-2xl border border-glass-border p-8 mb-20"
          >
            <h2 className="text-2xl font-bold text-foreground mb-8 text-center">
              Powered by Leading Cheminformatics Tools
            </h2>
            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-6">
              {INTEGRATIONS.map((integration, index) => {
                const Icon = integration.icon;
                return (
                  <motion.div
                    key={integration.name}
                    initial={{ opacity: 0, scale: 0.8 }}
                    animate={{ opacity: 1, scale: 1 }}
                    transition={{ duration: 0.5, delay: 0.5 + index * 0.1 }}
                    whileHover={{ scale: 1.05, y: -5 }}
                    className="flex flex-col items-center p-4 rounded-xl glass border border-glass-border hover:border-blue-500/50 transition-all cursor-pointer"
                  >
                    <div className={`w-12 h-12 rounded-xl bg-gradient-to-br ${integration.color} flex items-center justify-center mb-3 shadow-lg`}>
                      <Icon className="w-6 h-6 text-white" />
                    </div>
                    <p className="font-semibold text-sm text-foreground text-center">
                      {integration.name}
                    </p>
                    <p className="text-xs text-muted-foreground text-center mt-1">
                      {integration.description}
                    </p>
                  </motion.div>
                );
              })}
            </div>
          </motion.div>

          {/* Features Grid */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8 mb-20">
            {FEATURES.map((feature, index) => {
              const Icon = feature.icon;
              return (
                <motion.div
                  key={feature.title}
                  initial={{ opacity: 0, y: 40 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ duration: 0.8, delay: 0.6 + index * 0.1 }}
                  whileHover={{ y: -8 }}
                  className="glass-strong rounded-2xl p-8 border border-glass-border hover:border-blue-500/50 transition-all group"
                >
                  <div className={`w-14 h-14 rounded-xl bg-gradient-to-br ${feature.gradient} flex items-center justify-center mb-6 shadow-lg group-hover:scale-110 transition-transform`}>
                    <Icon className="w-7 h-7 text-white" />
                  </div>
                  <h3 className="text-xl font-semibold text-foreground mb-3">
                    {feature.title}
                  </h3>
                  <p className="text-muted-foreground leading-relaxed">
                    {feature.description}
                  </p>
                </motion.div>
              );
            })}
          </div>

          {/* CTA */}
          <motion.div
            initial={{ opacity: 0, y: 40 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.8, delay: 0.9 }}
            className="text-center glass-strong rounded-2xl border border-glass-border p-12"
          >
            <h2 className="text-3xl font-bold text-foreground mb-4">
              Ready to accelerate your drug discovery?
            </h2>
            <p className="text-muted-foreground mb-8 text-lg">
              Join researchers worldwide using PharmaSight for breakthrough discoveries
            </p>
            <Button
              size="lg"
              onClick={() => window.location.href = getLoginUrl()}
              className="bg-gradient-to-r from-blue-600 to-indigo-600 hover:from-blue-700 hover:to-indigo-700 text-lg px-10 py-6"
            >
              Sign In to PharmaSight
            </Button>
          </motion.div>
        </section>
      </div>
    );
  }

  // Authenticated view
  return (
    <div className="min-h-screen bg-background">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Chatbot */}
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            className="lg:col-span-2 h-[600px]"
          >
            <ChatbotInterface />
          </motion.div>

          {/* Quick Stats */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            className="space-y-4"
          >
            <div className="glass-strong rounded-xl border border-glass-border p-6">
              <h3 className="font-semibold text-foreground mb-4">Quick Access</h3>
              <div className="space-y-2">
                <Link href="/admin/dashboard">
                  <Button variant="outline" className="w-full justify-start glass border-glass-border">
                    View All Analogs
                  </Button>
                </Link>
                <Link href="/analytics">
                  <Button variant="outline" className="w-full justify-start glass border-glass-border">
                    Analytics Dashboard
                  </Button>
                </Link>
                <Link href="/testing">
                  <Button variant="outline" className="w-full justify-start glass border-glass-border">
                    Run Analysis
                  </Button>
                </Link>
              </div>
            </div>

            <div className="glass-strong rounded-xl border border-glass-border p-6 bg-gradient-to-br from-blue-500/10 to-indigo-500/10">
              <h3 className="font-semibold text-foreground mb-2">Latest Discoveries</h3>
              <p className="text-sm text-muted-foreground mb-4">
                Check the analytics dashboard for the latest analog discoveries and their properties.
              </p>
              <Link href="/analytics">
                <Button size="sm" className="w-full bg-gradient-to-r from-blue-600 to-indigo-600">
                  View Analytics
                </Button>
              </Link>
            </div>
          </motion.div>
        </div>
      </div>
    </div>
  );
}
