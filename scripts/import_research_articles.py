#!/usr/bin/env python3
"""
Import Research Articles into PostgreSQL Database
Reads from RESEARCH_ARTICLES_DATABASE.json and inserts into research_articles table
"""

import json
import psycopg2
from psycopg2.extras import execute_values
from datetime import datetime
import os

# Database connection from environment
DATABASE_URL = os.getenv('DATABASE_URL', 'postgresql://pharmasight_user:pharmasight_pass_2024@localhost:5432/pharmasight_db')

def import_research_articles():
    """Import all research articles from JSON to PostgreSQL"""

    # Load JSON data
    print("📂 Loading RESEARCH_ARTICLES_DATABASE.json...")
    with open('RESEARCH_ARTICLES_DATABASE.json', 'r') as f:
        articles_data = json.load(f)

    print(f"✅ Loaded {len(articles_data)} research articles")

    # Connect to database
    print(f"\n🔌 Connecting to PostgreSQL...")
    conn = psycopg2.connect(DATABASE_URL)
    cur = conn.cursor()

    # Prepare data for insertion
    insert_data = []
    skipped = 0

    for article in articles_data:
        try:
            # Extract fields
            pmid = article.get('pmid', article.get('PMID'))
            doi = article.get('doi', article.get('DOI'))
            title = article.get('title', '')
            abstract = article.get('abstract', '')

            # Authors (convert list to array)
            authors = article.get('authors', [])
            if isinstance(authors, str):
                authors = [a.strip() for a in authors.split(',')]

            journal = article.get('journal', '')

            # Publication date
            pub_date_str = article.get('publication_date', article.get('date'))
            pub_date = None
            year = article.get('year')

            if pub_date_str:
                try:
                    pub_date = datetime.strptime(pub_date_str, '%Y-%m-%d').date()
                    if not year:
                        year = pub_date.year
                except:
                    pass

            if not year and pub_date_str:
                try:
                    year = int(pub_date_str[:4])
                except:
                    pass

            # Keywords and MeSH terms
            keywords = article.get('keywords', [])
            if isinstance(keywords, str):
                keywords = [k.strip() for k in keywords.split(',')]

            mesh_terms = article.get('mesh_terms', article.get('MeSH_terms', []))
            if isinstance(mesh_terms, str):
                mesh_terms = [m.strip() for m in mesh_terms.split(',')]

            # Relevance scoring
            relevance_score = article.get('relevance_score', 0.75)
            cited_by_count = article.get('cited_by_count', article.get('citations', 0))

            insert_data.append((
                pmid, doi, title, abstract, authors, journal,
                pub_date, year, keywords, mesh_terms,
                relevance_score, cited_by_count
            ))

        except Exception as e:
            print(f"⚠️  Skipping article {article.get('pmid', 'unknown')}: {e}")
            skipped += 1
            continue

    # Bulk insert
    print(f"\n💾 Inserting {len(insert_data)} research articles into database...")

    insert_query = """
        INSERT INTO research_articles (
            pmid, doi, title, abstract, authors, journal,
            publication_date, year, keywords, mesh_terms,
            relevance_score, cited_by_count
        ) VALUES %s
        ON CONFLICT (pmid) DO UPDATE SET
            updated_at = CURRENT_TIMESTAMP,
            title = EXCLUDED.title,
            abstract = EXCLUDED.abstract,
            relevance_score = EXCLUDED.relevance_score
    """

    execute_values(cur, insert_query, insert_data)
    conn.commit()

    print(f"✅ Successfully imported {len(insert_data)} research articles")
    if skipped > 0:
        print(f"⚠️  Skipped {skipped} articles due to errors")

    # Query summary statistics
    cur.execute("SELECT COUNT(*) FROM research_articles")
    total_count = cur.fetchone()[0]

    cur.execute("SELECT COUNT(DISTINCT year) FROM research_articles WHERE year IS NOT NULL")
    year_span = cur.fetchone()[0]

    cur.execute("SELECT MIN(year), MAX(year) FROM research_articles WHERE year IS NOT NULL")
    year_range = cur.fetchone()

    cur.execute("SELECT COUNT(*) FROM research_articles WHERE doi IS NOT NULL")
    with_doi = cur.fetchone()[0]

    print(f"\n📊 Database Summary:")
    print(f"   Total articles in database: {total_count}")
    print(f"   Year range: {year_range[0]} - {year_range[1]}")
    print(f"   Unique years: {year_span}")
    print(f"   Articles with DOI: {with_doi} ({100*with_doi//total_count}%)")

    # Close connection
    cur.close()
    conn.close()

    print("\n🎉 Import complete!")

if __name__ == '__main__':
    import sys

    if not os.path.exists('RESEARCH_ARTICLES_DATABASE.json'):
        print("❌ Error: RESEARCH_ARTICLES_DATABASE.json not found")
        print("   Run this script from the pharmasight-platform root directory")
        sys.exit(1)

    try:
        import_research_articles()
    except Exception as e:
        print(f"\n❌ Import failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
