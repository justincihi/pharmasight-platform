# PharmaSight™ Research Article Database

## Overview

This database contains all research articles scanned by the PharmaSight™ autonomous research engine. It serves as a transparent record of the literature review process and supports intellectual property documentation.

## Database Statistics

- **Total Articles:** 26
- **Articles with DOI:** 25
- **Articles with PMID:** 22
- **Unique Journals:** 26
- **Year Range:** 1985 - 2024
- **Last Updated:** 2025-11-10T19:47:52.332774

## Source Databases

- PubMed

## Research Goals Covered

- Novel Anxiolytic Development from Natural Products
- Psychedelic Analog Development for Depression and Addiction
- ketamine NMDA antagonist rapid antidepressant
- muscimol GABA-A agonist sleep anxiety
- psilocybin 5-HT2A receptor depression treatment
- GABA-A Receptor Subtype-Selective Agonist Development
- MDMA PTSD therapy serotonin oxytocin
- NMDA Receptor Subtype-Selective Antagonist for Rapid Antidepressant Action
- Entactogen Development for PTSD Treatment
- Novel 5-HT2A Partial Agonist Discovery for Treatment-Resistant Depression

## File Formats

1. **JSON Database:** `RESEARCH_ARTICLES_DATABASE.json`
   - Complete article metadata including abstracts
   - Machine-readable format for programmatic access
   - Full search and filtering capabilities

2. **CSV Spreadsheet:** `RESEARCH_ARTICLES_DATABASE.csv`
   - Simple spreadsheet format for easy viewing
   - Compatible with Excel, Google Sheets, etc.
   - Includes: Article ID, Title, Authors, DOI, PMID, Journal, Year, URL, Keywords

## Usage

### Viewing the Database

**Option 1: CSV Spreadsheet (Easiest)**
- Open `RESEARCH_ARTICLES_DATABASE.csv` in Excel or Google Sheets
- Sort and filter by any column
- Search for specific articles or keywords

**Option 2: JSON Database (Programmatic)**
```python
import json

with open('RESEARCH_ARTICLES_DATABASE.json', 'r') as f:
    database = json.load(f)

# Access all articles
articles = database['articles']

# Search by DOI
doi_to_find = "10.1038/example"
article = next((a for a in articles if a.get('doi') == doi_to_find), None)
```

## Citation

When referencing this database, please use:

```
PharmaSight™ Research Article Database. Available at: 
https://github.com/justincihi/pharmasight-platform/blob/analog-discoveries-ip-protected/RESEARCH_ARTICLES_DATABASE.json
Last updated: 2025-11-10T19:47:52.332774
```

## Article Fields

Each article entry contains:

- **article_id:** Unique identifier (ART-00001, ART-00002, etc.)
- **title:** Full article title
- **authors:** List of authors
- **doi:** Digital Object Identifier
- **pmid:** PubMed ID
- **journal:** Journal name
- **year:** Publication year
- **abstract:** Article abstract (truncated to 500 chars)
- **url:** Direct link to article
- **keywords:** Relevant keywords
- **relevance_score:** Relevance to research goals (0-100)
- **scan_date:** Date article was scanned
- **source_database:** Database source (PubMed, arXiv, etc.)
- **research_goal:** Associated research goal

## Maintenance

This database is automatically updated by the PharmaSight™ autonomous research engine. Manual additions can be made using the ResearchArticleDatabase API.

## License

This database is maintained for research and intellectual property documentation purposes. Individual articles remain under their respective publishers' copyrights.

---

**Generated:** 2025-11-11 19:28:12 UTC
