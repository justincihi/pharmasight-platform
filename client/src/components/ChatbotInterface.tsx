import { useState, useRef, useEffect } from "react";
import { useLocation } from "wouter";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Badge } from "@/components/ui/badge";
import { toast } from "sonner";
import { Streamdown } from "streamdown";
import {
  MessageCircle,
  Send,
  Loader2,
  ExternalLink,
  ChevronRight,
  Database,
  Beaker,
  ShieldCheck,
  TrendingUp,
} from "lucide-react";

interface CitedCompound {
  id: number;
  compoundId: string;
  compoundName: string;
  confidenceScore: number;
  patentStatus: string;
  safetyScore: number;
  efficacyScore: number;
}

interface Message {
  id: string;
  role: "user" | "assistant";
  content: string;
  timestamp: Date;
  citedCompounds?: CitedCompound[];
  provider?: string;
}

const SAMPLE_QUESTIONS = [
  "Which analogs have the highest confidence scores?",
  "Show me all patent-free compounds with safety score above 80",
  "Compare the top 3 analogs by efficacy",
  "Which compounds have been flagged by ADMET screening?",
  "What are the best candidates for NMDA receptor docking?",
  "List all ketamine analogs discovered this month",
];

function PatentBadge({ status }: { status: string }) {
  const isPatentFree = status?.toLowerCase().includes("free") || status?.toLowerCase().includes("clear");
  return (
    <Badge
      variant="outline"
      className={`text-xs ${isPatentFree ? "border-emerald-500 text-emerald-600" : "border-amber-500 text-amber-600"}`}
    >
      {isPatentFree ? "Patent-free" : "Patented"}
    </Badge>
  );
}

function CompoundCitationCard({
  compound,
  onNavigate,
}: {
  compound: CitedCompound;
  onNavigate: (id: number) => void;
}) {
  return (
    <button
      onClick={() => onNavigate(compound.id)}
      className="w-full text-left group flex items-center gap-3 p-3 rounded-lg border border-border bg-card hover:bg-accent/50 hover:border-primary/40 transition-all duration-150"
    >
      <div className="flex-shrink-0 w-8 h-8 rounded-full bg-primary/10 flex items-center justify-center">
        <Beaker className="w-4 h-4 text-primary" />
      </div>
      <div className="flex-1 min-w-0">
        <div className="flex items-center gap-2 flex-wrap">
          <span className="font-medium text-sm text-foreground truncate">{compound.compoundName}</span>
          <PatentBadge status={compound.patentStatus} />
        </div>
        <div className="flex items-center gap-3 mt-1">
          <span className="text-xs text-muted-foreground flex items-center gap-1">
            <TrendingUp className="w-3 h-3" />
            {compound.confidenceScore}% confidence
          </span>
          <span className="text-xs text-muted-foreground flex items-center gap-1">
            <ShieldCheck className="w-3 h-3" />
            Safety {compound.safetyScore}/100
          </span>
        </div>
      </div>
      <ChevronRight className="w-4 h-4 text-muted-foreground group-hover:text-primary transition-colors flex-shrink-0" />
    </button>
  );
}

export function ChatbotInterface() {
  const [, setLocation] = useLocation();
  const [messages, setMessages] = useState<Message[]>([
    {
      id: "1",
      role: "assistant",
      content:
        "Hello! I'm **PharmaSight AI**, your drug discovery assistant. I have **live access** to your analog database, ADMET screening results, and docking data.\n\nAsk me anything — I can compare compounds, rank by safety/efficacy, explain SAR relationships, or surface patent-free candidates.",
      timestamp: new Date(),
    },
  ]);
  const [input, setInput] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const messagesEndRef = useRef<HTMLDivElement>(null);

  const chatMutation = trpc.chat.send.useMutation({
    onSuccess: (response: any) => {
      const assistantMessage: Message = {
        id: (Date.now() + 1).toString(),
        role: "assistant",
        content: response.content || "No response",
        timestamp: new Date(),
        citedCompounds: response.citedCompounds || [],
        provider: response.provider,
      };
      setMessages((prev) => [...prev, assistantMessage]);
      setIsLoading(false);
    },
    onError: (error: any) => {
      toast.error(`Chat error: ${error.message}`);
      setIsLoading(false);
    },
  });

  // Auto-scroll to bottom on new messages
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages, isLoading]);

  const handleSendMessage = async (text: string) => {
    if (!text.trim() || isLoading) return;

    const userMessage: Message = {
      id: Date.now().toString(),
      role: "user",
      content: text,
      timestamp: new Date(),
    };
    setMessages((prev) => [...prev, userMessage]);
    setInput("");
    setIsLoading(true);

    chatMutation.mutate({
      message: text,
      provider: "openai",
      history: messages
        .filter((m) => m.id !== "1")
        .map((m) => ({ role: m.role, content: m.content })),
    });
  };

  const handleNavigateToCompound = (id: number) => {
    setLocation(`/admin/analog/${id}`);
  };

  return (
    <div className="flex flex-col h-full bg-card rounded-xl border border-border overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-primary to-primary/80 text-primary-foreground px-4 py-3 flex items-center gap-3">
        <div className="w-8 h-8 rounded-full bg-white/20 flex items-center justify-center flex-shrink-0">
          <MessageCircle className="w-4 h-4" />
        </div>
        <div className="flex-1 min-w-0">
          <h2 className="font-semibold text-sm leading-tight">PharmaSight AI Assistant</h2>
          <p className="text-xs text-primary-foreground/70 leading-tight">Live database access · ADMET · Docking</p>
        </div>
        <div className="flex items-center gap-1 text-xs text-primary-foreground/70">
          <Database className="w-3 h-3" />
          <span>Live</span>
        </div>
      </div>

      {/* Messages Area */}
      <div className="flex-1 overflow-y-auto p-4 space-y-4 min-h-0">
        {messages.map((message) => (
          <div
            key={message.id}
            className={`flex flex-col gap-2 ${message.role === "user" ? "items-end" : "items-start"}`}
          >
            {/* Bubble */}
            <div
              className={`max-w-[85%] rounded-2xl px-4 py-3 ${
                message.role === "user"
                  ? "bg-primary text-primary-foreground rounded-br-sm"
                  : "bg-muted text-foreground rounded-bl-sm"
              }`}
            >
              {message.role === "assistant" ? (
                <div className="text-sm prose prose-sm dark:prose-invert max-w-none [&_table]:text-xs [&_th]:py-1 [&_td]:py-1">
                  <Streamdown>{message.content}</Streamdown>
                </div>
              ) : (
                <p className="text-sm">{message.content}</p>
              )}
              <p
                className={`text-xs mt-1.5 ${
                  message.role === "user" ? "text-primary-foreground/60" : "text-muted-foreground"
                }`}
              >
                {message.timestamp.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" })}
                {message.provider && message.role === "assistant" && (
                  <span className="ml-2 opacity-60">· {message.provider}</span>
                )}
              </p>
            </div>

            {/* Cited compound cards */}
            {message.role === "assistant" &&
              message.citedCompounds &&
              message.citedCompounds.length > 0 && (
                <div className="w-full max-w-[85%] space-y-1.5">
                  <p className="text-xs text-muted-foreground flex items-center gap-1 ml-1">
                    <ExternalLink className="w-3 h-3" />
                    Referenced compounds — click to view details
                  </p>
                  {message.citedCompounds.map((compound) => (
                    <CompoundCitationCard
                      key={compound.id}
                      compound={compound}
                      onNavigate={handleNavigateToCompound}
                    />
                  ))}
                </div>
              )}
          </div>
        ))}

        {/* Loading indicator */}
        {isLoading && (
          <div className="flex items-start gap-2">
            <div className="bg-muted rounded-2xl rounded-bl-sm px-4 py-3 flex items-center gap-2">
              <Loader2 className="w-4 h-4 animate-spin text-muted-foreground" />
              <span className="text-sm text-muted-foreground">Querying database…</span>
            </div>
          </div>
        )}
        <div ref={messagesEndRef} />
      </div>

      {/* Sample Questions — only show on first message */}
      {messages.length === 1 && (
        <div className="border-t border-border px-4 py-3 bg-muted/30">
          <p className="text-xs font-medium text-muted-foreground mb-2">Try asking:</p>
          <div className="grid grid-cols-1 gap-1.5">
            {SAMPLE_QUESTIONS.slice(0, 4).map((question, idx) => (
              <button
                key={idx}
                onClick={() => handleSendMessage(question)}
                className="text-left text-xs px-3 py-2 rounded-lg border border-border bg-background hover:bg-accent hover:border-primary/30 transition-all duration-150 text-foreground"
              >
                {question}
              </button>
            ))}
          </div>
        </div>
      )}

      {/* Input Area */}
      <div className="border-t border-border p-3 bg-background">
        <div className="flex gap-2">
          <Input
            placeholder="Ask about analogs, ADMET results, docking scores…"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter" && !e.shiftKey && !isLoading) {
                e.preventDefault();
                handleSendMessage(input);
              }
            }}
            disabled={isLoading}
            className="flex-1 text-sm"
          />
          <Button
            size="sm"
            onClick={() => handleSendMessage(input)}
            disabled={isLoading || !input.trim()}
            className="px-3"
          >
            {isLoading ? (
              <Loader2 className="w-4 h-4 animate-spin" />
            ) : (
              <Send className="w-4 h-4" />
            )}
          </Button>
        </div>
      </div>
    </div>
  );
}
