"""
Data Export Module for PharmaSight Platform
Supports CSV, Excel, and PDF export formats for research findings, compounds, and analytics
"""

import csv
import json
from datetime import datetime
from io import BytesIO, StringIO
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.lib.enums import TA_CENTER, TA_LEFT

class DataExporter:
    """Handle data export in multiple formats"""
    
    def __init__(self):
        self.styles = getSampleStyleSheet()
        self.title_style = ParagraphStyle(
            'CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1e40af'),
            spaceAfter=30,
            alignment=TA_CENTER
        )
        
    def export_to_csv(self, data, filename="export.csv"):
        """Export data to CSV format"""
        output = StringIO()
        
        if isinstance(data, list) and len(data) > 0:
            # Extract headers from first item
            headers = list(data[0].keys())
            writer = csv.DictWriter(output, fieldnames=headers)
            writer.writeheader()
            writer.writerows(data)
        elif isinstance(data, dict):
            # Single record
            writer = csv.DictWriter(output, fieldnames=data.keys())
            writer.writeheader()
            writer.writerow(data)
        
        return output.getvalue()
    
    def export_to_excel(self, data, filename="export.xlsx", sheet_name="Data"):
        """Export data to Excel format"""
        output = BytesIO()
        
        if isinstance(data, list):
            df = pd.DataFrame(data)
        elif isinstance(data, dict):
            df = pd.DataFrame([data])
        else:
            df = pd.DataFrame()
        
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            # Auto-adjust column widths
            worksheet = writer.sheets[sheet_name]
            for idx, col in enumerate(df.columns):
                max_length = max(
                    df[col].astype(str).apply(len).max(),
                    len(str(col))
                )
                worksheet.column_dimensions[chr(65 + idx)].width = min(max_length + 2, 50)
        
        output.seek(0)
        return output.getvalue()
    
    def export_compound_to_pdf(self, compound_data, filename="compound_report.pdf"):
        """Export compound analysis to PDF report"""
        output = BytesIO()
        doc = SimpleDocTemplate(output, pagesize=letter, topMargin=0.5*inch, bottomMargin=0.5*inch)
        story = []
        
        # Title
        title = Paragraph(f"<b>Compound Analysis Report</b>", self.title_style)
        story.append(title)
        story.append(Spacer(1, 0.3*inch))
        
        # Metadata
        metadata_style = self.styles['Normal']
        metadata = [
            f"<b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"<b>Platform:</b> PharmaSight™ v3.0.0",
            f"<b>Report Type:</b> Compound Analysis"
        ]
        for item in metadata:
            story.append(Paragraph(item, metadata_style))
        story.append(Spacer(1, 0.3*inch))
        
        # Compound Information
        story.append(Paragraph("<b>🧬 Compound Information</b>", self.styles['Heading2']))
        story.append(Spacer(1, 0.1*inch))
        
        compound_info = [
            ['Property', 'Value'],
            ['Compound Name', compound_data.get('name', 'N/A')],
            ['SMILES', compound_data.get('smiles', 'N/A')],
            ['Molecular Weight', f"{compound_data.get('molecular_weight', 'N/A')} g/mol"],
            ['LogP', str(compound_data.get('logP', 'N/A'))],
            ['TPSA', f"{compound_data.get('tpsa', 'N/A')} Ų"],
            ['Drug Likeness', f"{compound_data.get('drug_likeness', 'N/A')}%"],
        ]
        
        table = Table(compound_info, colWidths=[2.5*inch, 4*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1e40af')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
        ]))
        story.append(table)
        story.append(Spacer(1, 0.3*inch))
        
        # Therapeutic Information
        if compound_data.get('therapeutic_area'):
            story.append(Paragraph("<b>💊 Therapeutic Information</b>", self.styles['Heading2']))
            story.append(Spacer(1, 0.1*inch))
            
            therapeutic_info = [
                ['Property', 'Value'],
                ['Therapeutic Area', compound_data.get('therapeutic_area', 'N/A')],
                ['Development Status', compound_data.get('development_status', 'N/A')],
                ['Safety Score', f"{compound_data.get('safety_score', 'N/A')}%"],
                ['Efficacy Score', f"{compound_data.get('efficacy_score', 'N/A')}%"],
            ]
            
            table2 = Table(therapeutic_info, colWidths=[2.5*inch, 4*inch])
            table2.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#059669')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 12),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgreen),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
            ]))
            story.append(table2)
            story.append(Spacer(1, 0.3*inch))
        
        # Patent Information
        if compound_data.get('patent_status'):
            story.append(Paragraph("<b>📜 Patent Information</b>", self.styles['Heading2']))
            story.append(Spacer(1, 0.1*inch))
            
            patent_info = [
                ['Property', 'Value'],
                ['Patent Status', compound_data.get('patent_status', 'N/A')],
                ['Patent Number', compound_data.get('patent_number', 'N/A')],
                ['Expiration Date', compound_data.get('patent_expiration', 'N/A')],
            ]
            
            table3 = Table(patent_info, colWidths=[2.5*inch, 4*inch])
            table3.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#7c3aed')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 12),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lavender),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
            ]))
            story.append(table3)
        
        # Build PDF
        doc.build(story)
        output.seek(0)
        return output.getvalue()
    
    def export_research_findings_to_pdf(self, findings, filename="research_findings.pdf"):
        """Export research findings to PDF report"""
        output = BytesIO()
        doc = SimpleDocTemplate(output, pagesize=letter, topMargin=0.5*inch, bottomMargin=0.5*inch)
        story = []
        
        # Title
        title = Paragraph(f"<b>Research Findings Report</b>", self.title_style)
        story.append(title)
        story.append(Spacer(1, 0.2*inch))
        
        # Summary
        summary_style = self.styles['Normal']
        summary_text = f"""
        <b>Total Findings:</b> {len(findings)}<br/>
        <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br/>
        <b>Platform:</b> PharmaSight™ v3.0.0
        """
        story.append(Paragraph(summary_text, summary_style))
        story.append(Spacer(1, 0.3*inch))
        
        # Individual Findings
        for idx, finding in enumerate(findings, 1):
            story.append(Paragraph(f"<b>Finding #{idx}: {finding.get('title', 'Untitled')}</b>", 
                                 self.styles['Heading3']))
            story.append(Spacer(1, 0.1*inch))
            
            finding_data = [
                ['Property', 'Value'],
                ['Compound', finding.get('compound', 'N/A')],
                ['Therapeutic Area', finding.get('therapeutic_area', 'N/A')],
                ['Confidence', f"{finding.get('confidence', 0)}%"],
                ['Patent Potential', finding.get('patent_potential', 'N/A')],
                ['IP Status', finding.get('ip_status', 'N/A')],
                ['Estimated Value', finding.get('estimated_value', 'N/A')],
            ]
            
            table = Table(finding_data, colWidths=[2*inch, 4.5*inch])
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1e40af')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 10),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 8),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lightblue),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 1), (-1, -1), 9),
            ]))
            story.append(table)
            
            # Description
            if finding.get('description'):
                story.append(Spacer(1, 0.1*inch))
                desc_text = f"<b>Description:</b> {finding.get('description', '')}"
                story.append(Paragraph(desc_text, self.styles['Normal']))
            
            story.append(Spacer(1, 0.2*inch))
            
            # Page break after every 2 findings
            if idx % 2 == 0 and idx < len(findings):
                story.append(PageBreak())
        
        # Build PDF
        doc.build(story)
        output.seek(0)
        return output.getvalue()
    
    def export_analytics_to_pdf(self, analytics_data, filename="analytics_report.pdf"):
        """Export analytics dashboard to PDF report"""
        output = BytesIO()
        doc = SimpleDocTemplate(output, pagesize=letter, topMargin=0.5*inch, bottomMargin=0.5*inch)
        story = []
        
        # Title
        title = Paragraph(f"<b>Analytics Dashboard Report</b>", self.title_style)
        story.append(title)
        story.append(Spacer(1, 0.3*inch))
        
        # Research Productivity
        story.append(Paragraph("<b>📊 Research Productivity</b>", self.styles['Heading2']))
        story.append(Spacer(1, 0.1*inch))
        
        productivity_data = [
            ['Metric', 'Value'],
            ['Compounds Analyzed', str(analytics_data.get('compounds_analyzed', 0))],
            ['Analogs Generated', str(analytics_data.get('analogs_generated', 0))],
            ['Patents Identified', str(analytics_data.get('patents_identified', 0))],
        ]
        
        table1 = Table(productivity_data, colWidths=[3*inch, 3.5*inch])
        table1.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1e40af')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightblue),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
        ]))
        story.append(table1)
        story.append(Spacer(1, 0.3*inch))
        
        # Success Metrics
        story.append(Paragraph("<b>✨ Success Metrics</b>", self.styles['Heading2']))
        story.append(Spacer(1, 0.1*inch))
        
        success_data = [
            ['Metric', 'Value'],
            ['Hit Rate', f"{analytics_data.get('hit_rate', 0)}%"],
            ['Patent Success', f"{analytics_data.get('patent_success', 0)}%"],
            ['Time Saved', f"{analytics_data.get('time_saved', 0)} hours"],
        ]
        
        table2 = Table(success_data, colWidths=[3*inch, 3.5*inch])
        table2.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#059669')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightgreen),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
        ]))
        story.append(table2)
        
        # Build PDF
        doc.build(story)
        output.seek(0)
        return output.getvalue()


# Export functions for Flask integration
def export_data(data, format_type="csv", data_type="general", filename=None):
    """
    Main export function for Flask routes
    
    Args:
        data: Data to export (dict or list)
        format_type: 'csv', 'excel', or 'pdf'
        data_type: 'compound', 'research_findings', 'analytics', or 'general'
        filename: Optional custom filename
    
    Returns:
        Tuple of (data_bytes, mimetype, filename)
    """
    exporter = DataExporter()
    
    if format_type == "csv":
        output = exporter.export_to_csv(data, filename or "export.csv")
        return output, "text/csv", filename or "export.csv"
    
    elif format_type == "excel":
        output = exporter.export_to_excel(data, filename or "export.xlsx")
        return output, "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", filename or "export.xlsx"
    
    elif format_type == "pdf":
        if data_type == "compound":
            output = exporter.export_compound_to_pdf(data, filename or "compound_report.pdf")
        elif data_type == "research_findings":
            output = exporter.export_research_findings_to_pdf(data, filename or "research_findings.pdf")
        elif data_type == "analytics":
            output = exporter.export_analytics_to_pdf(data, filename or "analytics_report.pdf")
        else:
            # Generic PDF export (use compound format as default)
            output = exporter.export_compound_to_pdf(data, filename or "report.pdf")
        
        return output, "application/pdf", filename or "report.pdf"
    
    else:
        raise ValueError(f"Unsupported format type: {format_type}")

