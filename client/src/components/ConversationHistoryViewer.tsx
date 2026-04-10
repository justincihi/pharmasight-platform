import React, { useState, useEffect } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Badge } from '@/components/ui/badge';
import { Search, Download, Trash2, MessageSquare, Loader2 } from 'lucide-react';
import { toast } from 'sonner';

interface ConversationLog {
  id: string;
  date: Date;
  fileName: string;
  messageCount: number;
  topics: string[];
  duration: number; // in seconds
  preview: string;
}

export function ConversationHistoryViewer() {
  const [conversations, setConversations] = useState<ConversationLog[]>([]);
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedConversation, setSelectedConversation] = useState<ConversationLog | null>(null);
  const [conversationContent, setConversationContent] = useState<string>('');
  const [loading, setLoading] = useState(false);
  const [activeTab, setActiveTab] = useState('list');

  useEffect(() => {
    loadConversations();
  }, []);

  const loadConversations = async () => {
    setLoading(true);
    try {
      // In a real implementation, this would call a tRPC mutation to fetch logs
      // For now, we'll use mock data
      const mockConversations: ConversationLog[] = [
        {
          id: 'conv-1',
          date: new Date(Date.now() - 2 * 24 * 60 * 60 * 1000),
          fileName: 'conversation-admin-2026-04-08.md',
          messageCount: 12,
          topics: ['docking', 'admet', 'ketamine'],
          duration: 1245,
          preview: 'Discussion about ketamine analog docking and ADMET analysis...',
        },
        {
          id: 'conv-2',
          date: new Date(Date.now() - 5 * 24 * 60 * 60 * 1000),
          fileName: 'conversation-admin-2026-04-03.md',
          messageCount: 8,
          topics: ['discovery', 'scaffold', 'validation'],
          duration: 890,
          preview: 'Reviewed new discovery validation and scaffold library setup...',
        },
      ];
      setConversations(mockConversations);
    } catch (error) {
      toast.error('Failed to load conversations');
    } finally {
      setLoading(false);
    }
  };

  const handleViewConversation = async (conv: ConversationLog) => {
    setLoading(true);
    try {
      // In a real implementation, this would fetch the actual file content
      const mockContent = `# Conversation Log

**Date:** ${conv.date.toLocaleString()}
**Duration:** ${Math.round(conv.duration / 60)} minutes
**Total Messages:** ${conv.messageCount}
**Topics:** ${conv.topics.join(', ')}

---

## Message 1

**Role:** user
**Time:** ${conv.date.toLocaleString()}

**Content:**

Can you help me dock ketamine analogs against NMDA receptors?

---

## Message 2

**Role:** assistant
**Time:** ${new Date(conv.date.getTime() + 30000).toLocaleString()}

**Content:**

Of course! I can help you set up molecular docking simulations for ketamine analogs. Let me prepare the NMDA receptor structure and configure the docking parameters...

---`;

      setConversationContent(mockContent);
      setSelectedConversation(conv);
      setActiveTab('view');
    } catch (error) {
      toast.error('Failed to load conversation');
    } finally {
      setLoading(false);
    }
  };

  const handleDownloadConversation = async (conv: ConversationLog) => {
    try {
      // In a real implementation, this would trigger a download
      const element = document.createElement('a');
      element.setAttribute(
        'href',
        'data:text/markdown;charset=utf-8,' + encodeURIComponent(conversationContent)
      );
      element.setAttribute('download', conv.fileName);
      element.style.display = 'none';
      document.body.appendChild(element);
      element.click();
      document.body.removeChild(element);
      toast.success('Conversation downloaded');
    } catch (error) {
      toast.error('Failed to download conversation');
    }
  };

  const handleDeleteConversation = async (conv: ConversationLog) => {
    if (confirm('Are you sure you want to delete this conversation?')) {
      try {
        // In a real implementation, this would call a tRPC mutation
        setConversations(conversations.filter((c) => c.id !== conv.id));
        if (selectedConversation?.id === conv.id) {
          setSelectedConversation(null);
          setConversationContent('');
          setActiveTab('list');
        }
        toast.success('Conversation deleted');
      } catch (error) {
        toast.error('Failed to delete conversation');
      }
    }
  };

  const handleSearch = async () => {
    if (!searchQuery.trim()) {
      loadConversations();
      return;
    }

    setLoading(true);
    try {
      // In a real implementation, this would call a search tRPC mutation
      const filtered = conversations.filter(
        (conv) =>
          conv.preview.toLowerCase().includes(searchQuery.toLowerCase()) ||
          conv.topics.some((t) =>
            t.toLowerCase().includes(searchQuery.toLowerCase())
          )
      );
      setConversations(filtered);
    } catch (error) {
      toast.error('Search failed');
    } finally {
      setLoading(false);
    }
  };

  const filteredConversations = conversations.filter(
    (conv) =>
      conv.preview.toLowerCase().includes(searchQuery.toLowerCase()) ||
      conv.topics.some((t) =>
        t.toLowerCase().includes(searchQuery.toLowerCase())
      )
  );

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>Conversation History</CardTitle>
          <CardDescription>
            View and manage all chatbot conversations
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs value={activeTab} onValueChange={setActiveTab}>
            <TabsList className="grid w-full grid-cols-2">
              <TabsTrigger value="list">Conversations ({conversations.length})</TabsTrigger>
              <TabsTrigger value="view" disabled={!selectedConversation}>
                View Details
              </TabsTrigger>
            </TabsList>

            {/* List View */}
            <TabsContent value="list" className="space-y-4 mt-4">
              {/* Search */}
              <div className="flex gap-2">
                <div className="flex-1 relative">
                  <Search className="absolute left-3 top-3 h-4 w-4 text-muted-foreground" />
                  <Input
                    placeholder="Search by topic or content..."
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    className="pl-10"
                  />
                </div>
                <Button onClick={handleSearch} disabled={loading}>
                  {loading ? (
                    <Loader2 className="h-4 w-4 animate-spin" />
                  ) : (
                    'Search'
                  )}
                </Button>
              </div>

              {/* Conversation List */}
              {loading ? (
                <div className="flex items-center justify-center py-8">
                  <Loader2 className="h-6 w-6 animate-spin" />
                </div>
              ) : filteredConversations.length === 0 ? (
                <div className="text-center py-8 text-muted-foreground">
                  No conversations found
                </div>
              ) : (
                <div className="space-y-3">
                  {filteredConversations.map((conv) => (
                    <Card
                      key={conv.id}
                      className="cursor-pointer hover:bg-muted/50 transition-colors"
                      onClick={() => handleViewConversation(conv)}
                    >
                      <CardContent className="pt-6">
                        <div className="space-y-3">
                          {/* Header */}
                          <div className="flex items-start justify-between">
                            <div className="flex-1">
                              <div className="flex items-center gap-2">
                                <MessageSquare className="h-4 w-4 text-muted-foreground" />
                                <span className="font-semibold">
                                  {conv.date.toLocaleDateString()}
                                </span>
                              </div>
                              <p className="text-sm text-muted-foreground mt-1">
                                {conv.messageCount} messages • {Math.round(conv.duration / 60)} min
                              </p>
                            </div>
                            <div className="flex gap-2">
                              <Button
                                size="sm"
                                variant="ghost"
                                onClick={(e) => {
                                  e.stopPropagation();
                                  handleDownloadConversation(conv);
                                }}
                              >
                                <Download className="h-4 w-4" />
                              </Button>
                              <Button
                                size="sm"
                                variant="ghost"
                                onClick={(e) => {
                                  e.stopPropagation();
                                  handleDeleteConversation(conv);
                                }}
                              >
                                <Trash2 className="h-4 w-4" />
                              </Button>
                            </div>
                          </div>

                          {/* Topics */}
                          <div className="flex flex-wrap gap-2">
                            {conv.topics.map((topic) => (
                              <Badge key={topic} variant="secondary">
                                {topic}
                              </Badge>
                            ))}
                          </div>

                          {/* Preview */}
                          <p className="text-sm text-muted-foreground">
                            {conv.preview}
                          </p>
                        </div>
                      </CardContent>
                    </Card>
                  ))}
                </div>
              )}
            </TabsContent>

            {/* View Details */}
            <TabsContent value="view" className="space-y-4 mt-4">
              {selectedConversation && (
                <>
                  <div className="flex items-center justify-between">
                    <div>
                      <h3 className="font-semibold">{selectedConversation.fileName}</h3>
                      <p className="text-sm text-muted-foreground">
                        {selectedConversation.date.toLocaleString()}
                      </p>
                    </div>
                    <div className="flex gap-2">
                      <Button
                        onClick={() => handleDownloadConversation(selectedConversation)}
                      >
                        <Download className="mr-2 h-4 w-4" />
                        Download
                      </Button>
                      <Button
                        variant="destructive"
                        onClick={() => handleDeleteConversation(selectedConversation)}
                      >
                        <Trash2 className="mr-2 h-4 w-4" />
                        Delete
                      </Button>
                    </div>
                  </div>

                  {/* Content */}
                  <Card>
                    <CardContent className="pt-6">
                      <div className="bg-muted p-4 rounded-lg max-h-[600px] overflow-y-auto">
                        <pre className="text-sm whitespace-pre-wrap font-mono">
                          {conversationContent}
                        </pre>
                      </div>
                    </CardContent>
                  </Card>
                </>
              )}
            </TabsContent>
          </Tabs>
        </CardContent>
      </Card>
    </div>
  );
}
