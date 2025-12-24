import { Toaster } from "@/components/ui/sonner";
import { TooltipProvider } from "@/components/ui/tooltip";
import NotFound from "@/pages/NotFound";
import { Route, Switch } from "wouter";
import ErrorBoundary from "./components/ErrorBoundary";
import { ThemeProvider } from "./contexts/ThemeContext";
import { Navigation } from "./components/Navigation";
import Home from "./pages/Home";
import AdminDashboard from "./pages/AdminDashboard";
import Analytics from "./pages/Analytics";
import CompoundTesting from "./pages/CompoundTesting";
import BatchAnalysis from "./pages/BatchAnalysis";
import AnalogDetail from "./pages/AnalogDetail";
import SchedulerDashboard from "./pages/SchedulerDashboard";

function Router() {
  // make sure to consider if you need authentication for certain routes
  return (
    <>
      <Navigation />
      <Switch>
      <Route path={"/"} component={Home} />
      <Route path={"/admin/dashboard"} component={AdminDashboard} />
           <Route path="/analytics" component={Analytics} />
      <Route path="/testing" component={CompoundTesting} />
      <Route path="/batch" component={BatchAnalysis} />
      <Route path="/admin/analog/:id" component={AnalogDetail} />
      <Route path="/admin/scheduler" component={SchedulerDashboard} />
      <Route path={"/404"} component={NotFound} />
      {/* Final fallback route */}
      <Route component={NotFound} />
    </Switch>
    </>
  );
}

// NOTE: About Theme
// - First choose a default theme according to your design style (dark or light bg), than change color palette in index.css
//   to keep consistent foreground/background color across components
// - If you want to make theme switchable, pass `switchable` ThemeProvider and use `useTheme` hook

function App() {
  return (
    <ErrorBoundary>
      <ThemeProvider
        defaultTheme="light"
        // switchable
      >
        <TooltipProvider>
          <Toaster />
          <Router />
        </TooltipProvider>
      </ThemeProvider>
    </ErrorBoundary>
  );
}

export default App;
