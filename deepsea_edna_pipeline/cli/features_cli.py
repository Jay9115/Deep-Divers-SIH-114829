"""
Command-line interface for feature extraction pipeline.
"""

import click
from pathlib import Path
from ..module2_features.pipeline_features import FeaturePipeline


@click.group()
def features_cli():
    """Feature extraction pipeline commands."""
    pass


@features_cli.command()
@click.option('--input-fasta', '-i', required=True, type=click.Path(exists=True),
              help='Input FASTA file with ASV sequences')
@click.option('--output-dir', '-o', required=True, type=click.Path(),
              help='Output directory for features')
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Configuration YAML file')
@click.option('--run-id', help='Run identifier (default: auto-generated)')
def extract(input_fasta, output_dir, config, run_id):
    """Extract features from ASV sequences."""
    
    pipeline = FeaturePipeline(config_file=Path(config) if config else None)
    
    results = pipeline.run_feature_extraction(
        input_fasta=Path(input_fasta),
        output_dir=Path(output_dir),
        run_id=run_id
    )
    
    click.echo("Feature extraction completed successfully!")
    click.echo(f"K-mer features: {results['kmer_features']}")
    click.echo(f"FCGR features: {results['fcgr_features']}")
    click.echo(f"LSH index: {results['lsh_index']}")


if __name__ == '__main__':
    features_cli()