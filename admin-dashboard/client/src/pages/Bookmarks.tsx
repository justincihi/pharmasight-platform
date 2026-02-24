import { useState } from "react";
import { trpc } from "@/lib/trpc";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogFooter,
  DialogHeader,
  DialogTitle,
} from "@/components/ui/dialog";
import { Badge } from "@/components/ui/badge";
import { Bookmark, BookmarkX, Edit2, Trash2, ExternalLink, Filter, Search } from "lucide-react";
import { toast } from "sonner";
import { cn } from "@/lib/utils";
import { Link } from "wouter";
import { Navigation } from "@/components/Navigation";

interface BookmarkItem {
  id: number;
  userId: number;
  analogId: number | null;
  notificationId: number | null;
  title: string;
  notes: string | null;
  category: "high-priority" | "review-later" | "promising" | "archived";
  createdAt: Date | string;
  updatedAt: Date | string;
}

const categoryColors: Record<string, string> = {
  "high-priority": "bg-red-100 text-red-800 border-red-200",
  "review-later": "bg-blue-100 text-blue-800 border-blue-200",
  "promising": "bg-green-100 text-green-800 border-green-200",
  "archived": "bg-gray-100 text-gray-800 border-gray-200",
};

const categoryLabels: Record<string, string> = {
  "high-priority": "High Priority",
  "review-later": "Review Later",
  "promising": "Promising",
  "archived": "Archived",
};

export default function Bookmarks() {
  const [searchQuery, setSearchQuery] = useState("");
  const [filterCategory, setFilterCategory] = useState<string>("all");
  const [editingBookmark, setEditingBookmark] = useState<BookmarkItem | null>(null);
  const [editTitle, setEditTitle] = useState("");
  const [editNotes, setEditNotes] = useState("");
  const [editCategory, setEditCategory] = useState("");

  const utils = trpc.useUtils();

  // Fetch all bookmarks
  const { data: bookmarks = [], isLoading } = trpc.bookmarks.getAll.useQuery({ limit: 100 });

  // Update bookmark mutation
  const updateMutation = trpc.bookmarks.update.useMutation({
    onSuccess: () => {
      toast.success("Bookmark updated");
      utils.bookmarks.getAll.invalidate();
      setEditingBookmark(null);
    },
    onError: (error) => {
      toast.error("Failed to update bookmark: " + error.message);
    },
  });

  // Delete bookmark mutation
  const deleteMutation = trpc.bookmarks.delete.useMutation({
    onSuccess: () => {
      toast.success("Bookmark removed");
      utils.bookmarks.getAll.invalidate();
    },
    onError: (error) => {
      toast.error("Failed to delete bookmark: " + error.message);
    },
  });

  // Filter bookmarks
  const filteredBookmarks = (bookmarks as BookmarkItem[]).filter((bookmark) => {
    const matchesSearch = bookmark.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
      (bookmark.notes?.toLowerCase().includes(searchQuery.toLowerCase()) ?? false);
    const matchesCategory = filterCategory === "all" || bookmark.category === filterCategory;
    return matchesSearch && matchesCategory;
  });

  // Group bookmarks by category
  const groupedBookmarks = filteredBookmarks.reduce((acc, bookmark) => {
    const category = bookmark.category;
    if (!acc[category]) acc[category] = [];
    acc[category].push(bookmark);
    return acc;
  }, {} as Record<string, BookmarkItem[]>);

  const handleEdit = (bookmark: BookmarkItem) => {
    setEditingBookmark(bookmark);
    setEditTitle(bookmark.title);
    setEditNotes(bookmark.notes || "");
    setEditCategory(bookmark.category);
  };

  const handleSaveEdit = () => {
    if (!editingBookmark) return;
    
    updateMutation.mutate({
      bookmarkId: editingBookmark.id,
      title: editTitle,
      notes: editNotes,
      category: editCategory,
    });
  };

  const handleDelete = (bookmarkId: number) => {
    if (confirm("Are you sure you want to remove this bookmark?")) {
      deleteMutation.mutate({ bookmarkId });
    }
  };

  const formatDate = (dateInput: Date | string) => {
    const date = typeof dateInput === "string" ? new Date(dateInput) : dateInput;
    return date.toLocaleDateString("en-US", {
      month: "short",
      day: "numeric",
      year: "numeric",
    });
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <Navigation />
      
      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="flex items-center justify-between mb-8">
          <div>
            <h1 className="text-3xl font-bold text-gray-900 flex items-center gap-3">
              <Bookmark className="h-8 w-8 text-blue-600" />
              Saved Discoveries
            </h1>
            <p className="text-gray-600 mt-1">
              {bookmarks.length} saved item{bookmarks.length !== 1 ? "s" : ""}
            </p>
          </div>
        </div>

        {/* Filters */}
        <div className="flex flex-col sm:flex-row gap-4 mb-6">
          <div className="relative flex-1">
            <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-gray-400" />
            <Input
              placeholder="Search bookmarks..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              className="pl-10"
            />
          </div>
          <Select value={filterCategory} onValueChange={setFilterCategory}>
            <SelectTrigger className="w-full sm:w-48">
              <Filter className="h-4 w-4 mr-2" />
              <SelectValue placeholder="Filter by category" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="all">All Categories</SelectItem>
              <SelectItem value="high-priority">High Priority</SelectItem>
              <SelectItem value="promising">Promising</SelectItem>
              <SelectItem value="review-later">Review Later</SelectItem>
              <SelectItem value="archived">Archived</SelectItem>
            </SelectContent>
          </Select>
        </div>

        {/* Bookmarks List */}
        {isLoading ? (
          <div className="flex items-center justify-center py-12">
            <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600"></div>
          </div>
        ) : filteredBookmarks.length === 0 ? (
          <Card className="text-center py-12">
            <CardContent>
              <BookmarkX className="h-12 w-12 mx-auto text-gray-400 mb-4" />
              <h3 className="text-lg font-medium text-gray-900 mb-2">No bookmarks found</h3>
              <p className="text-gray-600 mb-4">
                {searchQuery || filterCategory !== "all"
                  ? "Try adjusting your search or filter"
                  : "Save important discoveries from notifications to see them here"}
              </p>
              <Link href="/">
                <Button variant="outline">
                  Go to Home
                </Button>
              </Link>
            </CardContent>
          </Card>
        ) : (
          <div className="space-y-6">
            {/* Show grouped by category */}
            {Object.entries(groupedBookmarks).map(([category, items]) => (
              <div key={category}>
                <h2 className="text-lg font-semibold text-gray-800 mb-3 flex items-center gap-2">
                  <Badge className={cn("font-normal", categoryColors[category])}>
                    {categoryLabels[category]}
                  </Badge>
                  <span className="text-sm text-gray-500">({items.length})</span>
                </h2>
                <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
                  {items.map((bookmark) => (
                    <Card key={bookmark.id} className="hover:shadow-md transition-shadow">
                      <CardHeader className="pb-2">
                        <div className="flex items-start justify-between">
                          <CardTitle className="text-base line-clamp-2">{bookmark.title}</CardTitle>
                          <Badge className={cn("ml-2 flex-shrink-0", categoryColors[bookmark.category])}>
                            {categoryLabels[bookmark.category]}
                          </Badge>
                        </div>
                        <CardDescription className="text-xs">
                          Saved on {formatDate(bookmark.createdAt)}
                        </CardDescription>
                      </CardHeader>
                      <CardContent>
                        {bookmark.notes && (
                          <p className="text-sm text-gray-600 mb-4 line-clamp-3">
                            {bookmark.notes}
                          </p>
                        )}
                        <div className="flex items-center gap-2">
                          {bookmark.analogId && (
                            <Link href={`/analog/${bookmark.analogId}`}>
                              <Button variant="outline" size="sm" className="flex-1">
                                <ExternalLink className="h-3 w-3 mr-1" />
                                View Details
                              </Button>
                            </Link>
                          )}
                          <Button
                            variant="ghost"
                            size="sm"
                            onClick={() => handleEdit(bookmark)}
                          >
                            <Edit2 className="h-3 w-3" />
                          </Button>
                          <Button
                            variant="ghost"
                            size="sm"
                            className="text-red-600 hover:text-red-700 hover:bg-red-50"
                            onClick={() => handleDelete(bookmark.id)}
                          >
                            <Trash2 className="h-3 w-3" />
                          </Button>
                        </div>
                      </CardContent>
                    </Card>
                  ))}
                </div>
              </div>
            ))}
          </div>
        )}

        {/* Edit Dialog */}
        <Dialog open={!!editingBookmark} onOpenChange={() => setEditingBookmark(null)}>
          <DialogContent>
            <DialogHeader>
              <DialogTitle>Edit Bookmark</DialogTitle>
              <DialogDescription>
                Update the details of your saved discovery
              </DialogDescription>
            </DialogHeader>
            <div className="space-y-4 py-4">
              <div>
                <label className="text-sm font-medium mb-1 block">Title</label>
                <Input
                  value={editTitle}
                  onChange={(e) => setEditTitle(e.target.value)}
                  placeholder="Bookmark title"
                />
              </div>
              <div>
                <label className="text-sm font-medium mb-1 block">Notes</label>
                <Textarea
                  value={editNotes}
                  onChange={(e) => setEditNotes(e.target.value)}
                  placeholder="Add your notes about this discovery..."
                  rows={4}
                />
              </div>
              <div>
                <label className="text-sm font-medium mb-1 block">Category</label>
                <Select value={editCategory} onValueChange={setEditCategory}>
                  <SelectTrigger>
                    <SelectValue />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="high-priority">High Priority</SelectItem>
                    <SelectItem value="promising">Promising</SelectItem>
                    <SelectItem value="review-later">Review Later</SelectItem>
                    <SelectItem value="archived">Archived</SelectItem>
                  </SelectContent>
                </Select>
              </div>
            </div>
            <DialogFooter>
              <Button variant="outline" onClick={() => setEditingBookmark(null)}>
                Cancel
              </Button>
              <Button onClick={handleSaveEdit} disabled={updateMutation.isPending}>
                {updateMutation.isPending ? "Saving..." : "Save Changes"}
              </Button>
            </DialogFooter>
          </DialogContent>
        </Dialog>
      </div>
    </div>
  );
}
