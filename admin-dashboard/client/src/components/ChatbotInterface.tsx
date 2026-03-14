import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Card } from "@/components/ui/card";
import { Send, MessageCircle, Loader2, ChevronDown } from "lucide-react";
import { trpc } from "@/lib/trpc";
import { toast } from "sonner";

interface Message {
  id: string;
  role: "user" | "assistant";
  content: string;
  timestamp: Date;
  provider?: string;
}

type LLMProvider = "openai" | "gemini" | "claude" | "perplexity" | "xai";

const PROVIDER_LABELS: Record<LLMProvider, string> = {
  openai: "OpenAI GPT-4o",
  gemini: "Google Gemini",
  claude: "Anthropic Claude",
  perplexity: "Perplexity Sonar",
  xai: "xAI Grok",
};

const SAMPLE_QUESTIONS = [
  "What are the top 5 analogs discovered this week?",
  "Show me patent-free compounds with >85% confidence",
  "Compare ADMET profiles of recent discoveries",
  "Which compounds have the highest market potential?",
  "Analyze safety scores for Ketamine analogs",
  "Generate a report on discovery trends",
];

export function ChatbotInterface() {
  const [messages, setMessages] = useState<Message[]>([
    {
      id: "1",
      role: "assistant",
      content: "Hello! I'm PharmaSight's AI assistant. I can help you explore analog discoveries, run cheminformatics analyses, and generate reports. What would you like to know?",
      timestamp: new Date(),
    },
  ]);
  const [input, setInput] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [selectedProvider, setSelectedProvider] = useState<LLMProvider>("openai");
  const [showProviderMenu, setShowProviderMenu] = useState(false);

  // Fetch which providers have API keys configured
  const { data: providersStatus } = trpc.chat.getProvidersStatus.useQuery(undefined, {
    retry: false,
  });

  const chatMutation = trpc.chat.send.useMutation({
    onSuccess: (response: any) => {
      const assistantMessage: Message = {
        id: (Date.now() + 1).toString(),
        role: "assistant",
        content: response.content || response.response || "No response",
        timestamp: new Date(),
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

  const handleSendMessage = async (text: string) => {
    if (!text.trim()) return;

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
      provider: selectedProvider,
      history: messages.map(m => ({ role: m.role, content: m.content })),
    });
  };

  const isProviderConfigured = (provider: LLMProvider): boolean => {
    if (!providersStatus) return provider === "openai";
    return providersStatus[provider]?.configured ?? false;
  };

  return (
    <div className="flex flex-col h-full bg-white rounded-lg border border-gray-200">
      {/* Header */}
      <div className="bg-gradient-to-r from-blue-600 to-blue-700 text-white p-4 rounded-t-lg">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <MessageCircle className="w-5 h-5" />
            <h2 className="font-semibold">PharmaSight AI Assistant</h2>
          </div>
          {/* Provider selector */}
          <div className="relative">
            <button
              className="flex items-center gap-1 text-xs bg-blue-500 hover:bg-blue-400 text-white px-2 py-1 rounded transition-colors"
              onClick={() => setShowProviderMenu((v) => !v)}
            >
              {PROVIDER_LABELS[selectedProvider]}
              <ChevronDown className="w-3 h-3" />
            </button>
            {showProviderMenu && (
              <div className="absolute right-0 mt-1 w-44 bg-white text-gray-800 border border-gray-200 rounded shadow-lg z-10">
                {(Object.keys(PROVIDER_LABELS) as LLMProvider[]).map((provider) => {
                  const configured = isProviderConfigured(provider);
                  return (
                    <button
                      key={provider}
                      className={`w-full text-left px-3 py-2 text-xs hover:bg-gray-100 flex items-center justify-between ${!configured ? "opacity-50 cursor-not-allowed" : ""}`}
                      onClick={() => {
                        if (!configured) {
                          toast.error(`${PROVIDER_LABELS[provider]} is not configured. Add the API key to your .env file.`);
                          return;
                        }
                        setSelectedProvider(provider);
                        setShowProviderMenu(false);
                      }}
                    >
                      <span>{PROVIDER_LABELS[provider]}</span>
                      {configured ? (
                        <span className="text-green-500 text-xs">✓</span>
                      ) : (
                        <span className="text-gray-400 text-xs">no key</span>
                      )}
                    </button>
                  );
                })}
              </div>
            )}
          </div>
        </div>
        <p className="text-blue-100 text-sm mt-1">
          Ask questions about your analog discoveries
        </p>
      </div>

      {/* Messages Area */}
      <div className="flex-1 overflow-y-auto p-4 space-y-4">
        {messages.map((message) => (
          <div
            key={message.id}
            className={`flex ${message.role === "user" ? "justify-end" : "justify-start"}`}
          >
            <div
              className={`max-w-xs lg:max-w-md px-4 py-2 rounded-lg ${
                message.role === "user"
                  ? "bg-blue-600 text-white rounded-br-none"
                  : "bg-gray-100 text-gray-900 rounded-bl-none"
              }`}
            >
              <p className="text-sm">{message.content}</p>
              <p
                className={`text-xs mt-1 ${
                  message.role === "user" ? "text-blue-100" : "text-gray-500"
                }`}
              >
                {message.timestamp.toLocaleTimeString([], {
                  hour: "2-digit",
                  minute: "2-digit",
                })}
                {message.provider && message.role === "assistant" && (
                  <span className="ml-1 opacity-70">· {PROVIDER_LABELS[message.provider as LLMProvider] ?? message.provider}</span>
                )}
              </p>
            </div>
          </div>
        ))}

        {isLoading && (
          <div className="flex justify-start">
            <div className="bg-gray-100 text-gray-900 px-4 py-2 rounded-lg rounded-bl-none">
              <Loader2 className="w-4 h-4 animate-spin" />
            </div>
          </div>
        )}
      </div>

      {/* Sample Questions */}
      {messages.length === 1 && (
        <div className="border-t border-gray-200 p-4 bg-gray-50">
          <p className="text-xs font-semibold text-gray-700 mb-3">Sample Questions</p>
          <div className="grid grid-cols-1 gap-2">
            {SAMPLE_QUESTIONS.slice(0, 3).map((question, idx) => (
              <Button
                key={idx}
                variant="outline"
                size="sm"
                className="justify-start text-left h-auto py-2 px-3"
                onClick={() => handleSendMessage(question)}
              >
                <span className="text-xs text-gray-700">{question}</span>
              </Button>
            ))}
          </div>
        </div>
      )}

      {/* Input Area */}
      <div className="border-t border-gray-200 p-4">
        <div className="flex gap-2">
          <Input
            placeholder="Ask about analogs, run tests, generate reports..."
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyPress={(e) => {
              if (e.key === "Enter" && !isLoading) {
                handleSendMessage(input);
              }
            }}
            disabled={isLoading}
            className="flex-1"
          />
          <Button
            size="sm"
            onClick={() => handleSendMessage(input)}
            disabled={isLoading || !input.trim()}
          >
            <Send className="w-4 h-4" />
          </Button>
        </div>
      </div>
    </div>
  );
}
