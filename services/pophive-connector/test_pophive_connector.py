"""
PopHIVE Connector Test Suite and Visualization Demo
====================================================

This script demonstrates the PopHIVE connector capabilities and generates
visualizations useful for the PharmaSight drug discovery platform.
"""

import sys
sys.path.insert(0, '/home/ubuntu/pharmasight-platform/services/pophive-connector')

from pophive_connector import PopHIVEConnector, get_diabetes_data, get_rsv_data
import pandas as pd
import matplotlib.pyplot as plt
import json
from datetime import datetime
import os

# Create output directory
OUTPUT_DIR = "/home/ubuntu/pharmasight-platform/data/pophive_exports"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def test_chronic_disease_data():
    """Test chronic disease data retrieval."""
    print("\n" + "=" * 70)
    print("TEST 1: Chronic Disease Data Retrieval")
    print("=" * 70)
    
    connector = PopHIVEConnector()
    
    # Get diabetes data
    diabetes_df = connector.get_chronic_disease_data(condition="diabetes")
    
    print(f"\n✓ Successfully retrieved diabetes data for {len(diabetes_df)} states")
    print(f"✓ Data year: {diabetes_df['year'].iloc[0]}")
    print(f"✓ Data source: {diabetes_df['source'].iloc[0]}")
    
    # Summary statistics
    print(f"\n📊 Summary Statistics:")
    print(f"   - Mean prevalence: {diabetes_df['value'].mean():.2f}%")
    print(f"   - Median prevalence: {diabetes_df['value'].median():.2f}%")
    print(f"   - Std deviation: {diabetes_df['value'].std():.2f}%")
    print(f"   - Min prevalence: {diabetes_df['value'].min():.2f}% ({diabetes_df.loc[diabetes_df['value'].idxmin(), 'geography']})")
    print(f"   - Max prevalence: {diabetes_df['value'].max():.2f}% ({diabetes_df.loc[diabetes_df['value'].idxmax(), 'geography']})")
    print(f"   - Total sample size: {diabetes_df['sample_size'].sum():,}")
    
    # Export to CSV
    csv_path = f"{OUTPUT_DIR}/diabetes_prevalence_by_state.csv"
    diabetes_df.to_csv(csv_path, index=False)
    print(f"\n💾 Data exported to: {csv_path}")
    
    return diabetes_df


def test_respiratory_disease_data():
    """Test respiratory disease data retrieval."""
    print("\n" + "=" * 70)
    print("TEST 2: Respiratory Disease (RSV) Data Retrieval")
    print("=" * 70)
    
    connector = PopHIVEConnector()
    
    # Get RSV data
    rsv_df = connector.get_respiratory_disease_data(disease="rsv")
    
    print(f"\n✓ Successfully retrieved RSV data for {len(rsv_df)} states")
    print(f"✓ Data date: {rsv_df['date'].iloc[0].strftime('%Y-%m-%d')}")
    print(f"✓ Data source: {rsv_df['source'].iloc[0]}")
    
    # Summary statistics
    print(f"\n📊 Summary Statistics:")
    print(f"   - Mean ED visit rate: {rsv_df['value'].mean():.2f}%")
    print(f"   - Max ED visit rate: {rsv_df['value'].max():.2f}% ({rsv_df.loc[rsv_df['value'].idxmax(), 'geography']})")
    print(f"   - Min ED visit rate: {rsv_df['value'].min():.2f}% ({rsv_df.loc[rsv_df['value'].idxmin(), 'geography']})")
    
    # Export to CSV
    csv_path = f"{OUTPUT_DIR}/rsv_surveillance_by_state.csv"
    rsv_df.to_csv(csv_path, index=False)
    print(f"\n💾 Data exported to: {csv_path}")
    
    return rsv_df


def test_diabetes_trends():
    """Test diabetes trends over time."""
    print("\n" + "=" * 70)
    print("TEST 3: Diabetes Prevalence Trends (2018-2025)")
    print("=" * 70)
    
    connector = PopHIVEConnector()
    
    # Get trends data
    trends_df = connector.get_diabetes_trends()
    
    print(f"\n✓ Successfully retrieved trend data for {len(trends_df)} records")
    
    # Separate by source
    hba1c_trends = trends_df[trends_df['source'] == 'Epic Cosmos: HbA1c']
    icd10_trends = trends_df[trends_df['source'] == 'Epic Cosmos: ICD10']
    
    print(f"\n📊 HbA1c-based Diabetes Trends:")
    for _, row in hba1c_trends.iterrows():
        print(f"   {row['year']}: {row['value']:.2f}% (n={row['sample_size']:,})")
    
    print(f"\n📊 ICD10-based Diabetes Trends:")
    for _, row in icd10_trends.iterrows():
        print(f"   {row['year']}: {row['value']:.2f}% (n={row['sample_size']:,})")
    
    # Export to CSV
    csv_path = f"{OUTPUT_DIR}/diabetes_trends_national.csv"
    trends_df.to_csv(csv_path, index=False)
    print(f"\n💾 Data exported to: {csv_path}")
    
    return trends_df


def test_high_prevalence_analysis():
    """Test high prevalence state identification."""
    print("\n" + "=" * 70)
    print("TEST 4: High Prevalence State Analysis")
    print("=" * 70)
    
    connector = PopHIVEConnector()
    
    # Get high prevalence states
    high_prev_df = connector.get_high_prevalence_states(threshold=8.0, top_n=15)
    
    print(f"\n✓ Identified {len(high_prev_df)} states with diabetes prevalence ≥8%")
    
    print(f"\n📊 High Prevalence States (Drug Intervention Priority):")
    print("-" * 50)
    for i, (_, row) in enumerate(high_prev_df.iterrows(), 1):
        print(f"   {i:2}. {row['geography']:20} {row['value']:.2f}% (n={row['sample_size']:,})")
    
    # Export to CSV
    csv_path = f"{OUTPUT_DIR}/high_prevalence_states.csv"
    high_prev_df.to_csv(csv_path, index=False)
    print(f"\n💾 Data exported to: {csv_path}")
    
    return high_prev_df


def test_pharmasight_integration():
    """Test full PharmaSight integration data package."""
    print("\n" + "=" * 70)
    print("TEST 5: PharmaSight Integration Data Package")
    print("=" * 70)
    
    connector = PopHIVEConnector()
    
    # Get full integration data
    integration_data = connector.get_pharmasight_integration_data()
    
    print(f"\n✓ Successfully generated integration data package")
    
    print(f"\n📊 Metadata:")
    print(f"   - Source: {integration_data['metadata']['source']}")
    print(f"   - Fetch timestamp: {integration_data['metadata']['fetch_timestamp']}")
    print(f"   - Data sources: {len(integration_data['metadata']['data_sources'])}")
    
    print(f"\n📊 Summary Statistics:")
    stats = integration_data['summary_statistics']
    print(f"   - States covered: {stats['total_states_covered']}")
    print(f"   - Avg diabetes prevalence: {stats['avg_diabetes_prevalence']}%")
    print(f"   - Max diabetes prevalence: {stats['max_diabetes_prevalence']}%")
    print(f"   - Min diabetes prevalence: {stats['min_diabetes_prevalence']}%")
    print(f"   - Total sample size: {stats['total_sample_size']:,}")
    
    # Export to JSON
    json_path = f"{OUTPUT_DIR}/pharmasight_integration_data.json"
    with open(json_path, 'w') as f:
        json.dump(integration_data, f, indent=2, default=str)
    print(f"\n💾 Data exported to: {json_path}")
    
    return integration_data


def create_visualizations(diabetes_df, rsv_df, trends_df, high_prev_df):
    """Create visualizations for PharmaSight dashboard."""
    print("\n" + "=" * 70)
    print("GENERATING VISUALIZATIONS")
    print("=" * 70)
    
    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # 1. Diabetes Prevalence by State (Top 20)
    fig, ax = plt.subplots(figsize=(12, 8))
    top_20 = diabetes_df.nlargest(20, 'value')
    colors = plt.cm.Reds(top_20['value'] / top_20['value'].max())
    bars = ax.barh(top_20['geography'], top_20['value'], color=colors)
    ax.set_xlabel('Diabetes Prevalence (%)', fontsize=12)
    ax.set_ylabel('State', fontsize=12)
    ax.set_title('Diabetes Prevalence by State (Top 20)\nData Source: PopHIVE / Epic Cosmos', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    
    # Add value labels
    for bar, val in zip(bars, top_20['value']):
        ax.text(val + 0.1, bar.get_y() + bar.get_height()/2, f'{val:.1f}%', 
                va='center', fontsize=9)
    
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/diabetes_prevalence_by_state.png", dpi=150, bbox_inches='tight')
    print(f"✓ Saved: diabetes_prevalence_by_state.png")
    plt.close()
    
    # 2. RSV ED Visits by State (Top 20)
    fig, ax = plt.subplots(figsize=(12, 8))
    top_20_rsv = rsv_df.nlargest(20, 'value')
    colors = plt.cm.Blues(top_20_rsv['value'] / top_20_rsv['value'].max())
    bars = ax.barh(top_20_rsv['geography'], top_20_rsv['value'], color=colors)
    ax.set_xlabel('RSV ED Visit Rate (%)', fontsize=12)
    ax.set_ylabel('State', fontsize=12)
    ax.set_title('RSV Emergency Department Visits by State (Top 20)\nData Source: PopHIVE / CDC NSSP', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    
    for bar, val in zip(bars, top_20_rsv['value']):
        ax.text(val + 0.01, bar.get_y() + bar.get_height()/2, f'{val:.2f}%', 
                va='center', fontsize=9)
    
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/rsv_ed_visits_by_state.png", dpi=150, bbox_inches='tight')
    print(f"✓ Saved: rsv_ed_visits_by_state.png")
    plt.close()
    
    # 3. Diabetes Trends Over Time
    fig, ax = plt.subplots(figsize=(10, 6))
    
    hba1c_trends = trends_df[trends_df['source'] == 'Epic Cosmos: HbA1c']
    icd10_trends = trends_df[trends_df['source'] == 'Epic Cosmos: ICD10']
    
    ax.plot(hba1c_trends['year'], hba1c_trends['value'], 'o-', linewidth=2, 
            markersize=8, label='HbA1c-based (Uncontrolled)', color='#e74c3c')
    ax.plot(icd10_trends['year'], icd10_trends['value'], 's-', linewidth=2, 
            markersize=8, label='ICD10-based (Diagnosed)', color='#3498db')
    
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Diabetes Prevalence (%)', fontsize=12)
    ax.set_title('National Diabetes Prevalence Trends (2018-2025)\nData Source: PopHIVE / Epic Cosmos', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/diabetes_trends_national.png", dpi=150, bbox_inches='tight')
    print(f"✓ Saved: diabetes_trends_national.png")
    plt.close()
    
    # 4. Combined Dashboard View
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Panel A: Diabetes Distribution
    ax1 = axes[0, 0]
    ax1.hist(diabetes_df['value'], bins=15, color='#e74c3c', edgecolor='white', alpha=0.8)
    ax1.axvline(diabetes_df['value'].mean(), color='black', linestyle='--', linewidth=2, label=f'Mean: {diabetes_df["value"].mean():.1f}%')
    ax1.set_xlabel('Diabetes Prevalence (%)', fontsize=11)
    ax1.set_ylabel('Number of States', fontsize=11)
    ax1.set_title('A. Distribution of Diabetes Prevalence', fontsize=12, fontweight='bold')
    ax1.legend()
    
    # Panel B: RSV Distribution
    ax2 = axes[0, 1]
    ax2.hist(rsv_df['value'], bins=15, color='#3498db', edgecolor='white', alpha=0.8)
    ax2.axvline(rsv_df['value'].mean(), color='black', linestyle='--', linewidth=2, label=f'Mean: {rsv_df["value"].mean():.2f}%')
    ax2.set_xlabel('RSV ED Visit Rate (%)', fontsize=11)
    ax2.set_ylabel('Number of States', fontsize=11)
    ax2.set_title('B. Distribution of RSV ED Visits', fontsize=12, fontweight='bold')
    ax2.legend()
    
    # Panel C: Top 10 High Prevalence States
    ax3 = axes[1, 0]
    top_10 = high_prev_df.head(10)
    colors = plt.cm.Reds(top_10['value'] / top_10['value'].max())
    bars = ax3.barh(top_10['geography'], top_10['value'], color=colors)
    ax3.set_xlabel('Diabetes Prevalence (%)', fontsize=11)
    ax3.set_title('C. Top 10 High Prevalence States', fontsize=12, fontweight='bold')
    ax3.invert_yaxis()
    
    # Panel D: Trends
    ax4 = axes[1, 1]
    hba1c_trends = trends_df[trends_df['source'] == 'Epic Cosmos: HbA1c']
    ax4.plot(hba1c_trends['year'], hba1c_trends['value'], 'o-', linewidth=2, 
            markersize=8, color='#e74c3c')
    ax4.fill_between(hba1c_trends['year'], hba1c_trends['value'], alpha=0.3, color='#e74c3c')
    ax4.set_xlabel('Year', fontsize=11)
    ax4.set_ylabel('Prevalence (%)', fontsize=11)
    ax4.set_title('D. National Diabetes Trend', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle('PopHIVE Data Dashboard for PharmaSight\nPopulation Health Insights for Drug Discovery', 
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/pophive_dashboard_overview.png", dpi=150, bbox_inches='tight')
    print(f"✓ Saved: pophive_dashboard_overview.png")
    plt.close()
    
    print(f"\n✓ All visualizations saved to: {OUTPUT_DIR}/")


def generate_pharmasight_use_cases():
    """Generate use case documentation for PharmaSight integration."""
    print("\n" + "=" * 70)
    print("PHARMASIGHT USE CASES")
    print("=" * 70)
    
    use_cases = """
# PopHIVE Data Use Cases for PharmaSight

## 1. Drug Safety Monitoring
- **Use Case**: Monitor chronic disease prevalence to identify potential drug safety signals
- **Data**: Diabetes/obesity prevalence by state over time
- **Application**: If a new diabetes drug is launched, track regional prevalence changes
- **Benefit**: Early detection of unexpected efficacy or safety patterns

## 2. Population PK/PD Modeling
- **Use Case**: Incorporate population health data into pharmacokinetic models
- **Data**: Disease prevalence, demographic distributions, comorbidity rates
- **Application**: Adjust dosing recommendations based on regional health profiles
- **Benefit**: More accurate population-level drug response predictions

## 3. Drug Interaction Analysis
- **Use Case**: Correlate medication use patterns with health outcomes
- **Data**: RSV/respiratory disease surveillance + chronic disease data
- **Application**: Identify potential drug-disease interactions in real-world populations
- **Benefit**: Enhanced DDI prediction with population context

## 4. Clinical Trial Site Selection
- **Use Case**: Identify optimal locations for clinical trials
- **Data**: High prevalence states for target conditions
- **Application**: Select trial sites with sufficient patient populations
- **Benefit**: Faster patient recruitment, more representative samples

## 5. Market Analysis for Drug Development
- **Use Case**: Assess market potential for new therapeutics
- **Data**: Disease prevalence trends, sample sizes by region
- **Application**: Estimate addressable patient populations
- **Benefit**: Data-driven drug development prioritization

## 6. Real-World Evidence Generation
- **Use Case**: Generate RWE for regulatory submissions
- **Data**: Multi-source health data (Epic Cosmos, BRFSS, Medicare)
- **Application**: Support drug efficacy claims with population-level data
- **Benefit**: Stronger regulatory evidence packages

## Data Sources Available
1. **Epic Cosmos** - Clinical EHR data (measured values)
2. **CDC BRFSS** - Survey-based health data
3. **Medicare FFS** - Claims-based data
4. **CDC NSSP** - Syndromic surveillance
5. **RESP-NET** - Respiratory disease monitoring
6. **National Wastewater Surveillance** - Environmental monitoring
7. **Google Trends** - Health search trends
8. **National Immunization Survey** - Vaccination data
"""
    
    # Save use cases document
    doc_path = f"{OUTPUT_DIR}/pharmasight_use_cases.md"
    with open(doc_path, 'w') as f:
        f.write(use_cases)
    
    print(use_cases)
    print(f"\n💾 Use cases document saved to: {doc_path}")


def main():
    """Run all tests and generate outputs."""
    print("\n" + "=" * 70)
    print("PopHIVE CONNECTOR TEST SUITE FOR PHARMASIGHT")
    print("=" * 70)
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output Directory: {OUTPUT_DIR}")
    
    # Run tests
    diabetes_df = test_chronic_disease_data()
    rsv_df = test_respiratory_disease_data()
    trends_df = test_diabetes_trends()
    high_prev_df = test_high_prevalence_analysis()
    integration_data = test_pharmasight_integration()
    
    # Create visualizations
    create_visualizations(diabetes_df, rsv_df, trends_df, high_prev_df)
    
    # Generate use cases
    generate_pharmasight_use_cases()
    
    # Final summary
    print("\n" + "=" * 70)
    print("TEST SUITE COMPLETE")
    print("=" * 70)
    print(f"\n📁 All outputs saved to: {OUTPUT_DIR}/")
    print("\nGenerated Files:")
    for f in os.listdir(OUTPUT_DIR):
        filepath = os.path.join(OUTPUT_DIR, f)
        size = os.path.getsize(filepath)
        print(f"   - {f} ({size:,} bytes)")
    
    print("\n✅ PopHIVE Connector is ready for PharmaSight integration!")


if __name__ == "__main__":
    main()
