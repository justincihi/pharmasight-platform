import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Badge } from '@/components/ui/badge';
import { Search, Download, Trash2, MessageSquare, Loader2 } from 'lucide-react';
import { toast } from 'sonner';
import { trpc } from '@/lib/trpc';

export function ConversationHistoryViewer() {
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedConversation, setSelectedConversation] = useState<string | null>(null);
  const utils = trpc.useUtils();

  // List conversations
  const { data: conversationsList, isLoading: isLoadingList } = 
    trpc.conversationLogger.listConversations.useQuery({ limit: 50, offset: 0 });

  // Get specific conversation
  const { data: selectedConvData, isLoading: isLoadingContent } = 
    trpc.conversationLogger.getConversation.useQuery(
      { filename: selectedConversation || "" },
      { enabled: !!selectedConversation }
    );

  // Search conversations
  const { data: searchResults, isLoading: isSearching } = 
    trpc.conversationLogger.searchConversations.useQuery(
      { query: searchQuery, limit: 20 },
      { enabled: searchQuery.length > 2 }
    );

  // Delete conversation
  const deleteConvMutation = trpc.conversationLogger.deleteConversation.useMutation({
    onSuccess: () => {
      toast.success("Conversation deleted");
      setSelectedConversation(null);
      utils.conversationLogger.listConversations.invalidate();
    },
    onError: (error: any) => {
      toast.error(`Delete failed: ${error.message}`);
    },
  });

  const handleExport = async () => {
    if (!selectedConversation) return;
    try {
      // Fetch export data
      const result = await utils.conversationLogger.exportConversation.fetch({ filename: selectedConversation });
      if (result.success && result.content && result.downloadName) {
        // Create download link
        const element = document.createElement("a");
        element.setAttribute("href", "data:text/markdown;charset=utf-8," + encodeURIComponent(result.content as string));
        element.setAttribute("download", result.downloadName as string);
        element.style.display = "none";
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
        toast.success("Conversation exported");
      }
    } catch (error: any) {
      toast.error(`Export failed: ${error.message}`);
    }
  };

  const handleSearch = (e: React.ChangeEvent<HTMLInputElement>) => {
    setSearchQuery(e.target.value);
  };

  const displayConversations = searchQuery.length > 2 && searchResults?.results 
    ? searchResults.results 
    : conversationsList?.conversations || [];

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>Conversation History</CardTitle>
          <CardDescription>
            View and manage all chatbot conversations
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Search Bar */}
          <div className="flex gap-2">
            <div className="flex-1 relative">
              <Search className="absolute left-3 top-3 h-4 w-4 text-muted-foreground" />
              <Input
                placeholder="Search conversations by topic or content..."
                value={searchQuery}
                onChange={handleSearch}
                className="pl-10"
              />
            </div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
            {/* Conversations List */}
            <div className="lg:col-span-1">
              <Card>
                <CardHeader>
                  <CardTitle className="text-sm">
                    {searchQuery.length > 2 ? "Search Results" : "Recent Conversations"}
                  </CardTitle>
                  <CardDescription className="text-xs">
                    {displayConversations.length} found
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="space-y-2 max-h-96 overflow-y-auto">
                    {isLoadingList || isSearching ? (
                      <div className="flex items-center justify-center py-4">
                        <Loader2 className="h-4 w-4 animate-spin" />
                      </div>
                    ) : !displayConversations || displayConversations.length === 0 ? (
                      <p className="text-sm text-muted-foreground">No conversations found</p>
                    ) : (
                      (displayConversations as any[]).map((conv: any) => (
                        <button
                          key={conv.filename}
                          onClick={() => setSelectedConversation(conv.filename)}
                          className={`w-full text-left p-2 rounded-md text-sm transition-colors ${
                            selectedConversation === conv.filename
                              ? "bg-primary text-primary-foreground"
                              : "hover:bg-muted"
                          }`}
                        >
                          <div className="truncate font-medium">{conv.filename}</div>
                          {conv.modifiedAt && (
                            <div className="text-xs opacity-70">
                              {new Date(conv.modifiedAt).toLocaleDateString()}
                            </div>
                          )}
                        </button>
                      ))
                    )}
                  </div>
                </CardContent>
              </Card>
            </div>

            {/* Conversation Content */}
            <div className="lg:col-span-2">
              {selectedConversation ? (
                <Card>
                  <CardHeader>
                    <div className="flex items-start justify-between">
                      <div>
                        <CardTitle className="text-base">{selectedConversation}</CardTitle>
                        <CardDescription>
                          Conversation details and messages
                        </CardDescription>
                      </div>
                      <div className="flex gap-2">
                        <Button
                          variant="outline"
                          size="sm"
                          onClick={handleExport}
                        >
                          <Download className="h-3 w-3 mr-1" />
                          Export
                        </Button>
                        <Button
                          variant="destructive"
                          size="sm"
                          onClick={() => deleteConvMutation.mutate({ filename: selectedConversation })}
                          disabled={deleteConvMutation.isPending}
                        >
                          <Trash2 className="h-3 w-3 mr-1" />
                          Delete
                        </Button>
                      </div>
                    </div>
                  </CardHeader>
                  <CardContent>
                    {isLoadingContent ? (
                      <div className="flex items-center justify-center py-8">
                        <Loader2 className="h-6 w-6 animate-spin" />
                      </div>
                    ) : selectedConvData?.success ? (
                      <div className="bg-muted p-4 rounded-lg max-h-[600px] overflow-y-auto">
                        <pre className="text-sm whitespace-pre-wrap font-mono">
                          {selectedConvData.content}
                        </pre>
                      </div>
                    ) : (
                      <p className="text-sm text-destructive">
                        {selectedConvData?.error || "Failed to load conversation"}
                      </p>
                    )}
                  </CardContent>
                </Card>
              ) : (
                <Card>
                  <CardContent className="flex flex-col items-center justify-center py-12">
                    <MessageSquare className="h-12 w-12 text-muted-foreground mb-4" />
                    <p className="text-muted-foreground">Select a conversation to view details</p>
                  </CardContent>
                </Card>
              )}
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
