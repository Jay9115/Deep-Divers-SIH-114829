"""
Command-line interface for ASV generation pipeline.
"""

import click
from pathlib import Path
from ..module1_qc_asv.pipeline_asv import ASVPipeline


@click.group()
def asv_cli():
    """ASV generation pipeline commands."""
    pass


@asv_cli.command()
@click.option('--input-dir', '-i', required=True, type=click.Path(exists=True),
              help='Input directory containing FASTQ files')
@click.option('--output-dir', '-o', required=True, type=click.Path(),
              help='Output directory for results')
@click.option('--metadata', '-m', type=click.Path(exists=True),
              help='Metadata file (TSV format)')
@click.option('--reference-db', '-r', type=click.Path(exists=True),
              help='Reference database for chimera detection')
@click.option('--negative-controls', '-n', multiple=True,
              help='Negative control sample names')
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Configuration YAML file')
@click.option('--run-id', help='Run identifier (default: auto-generated)')
def run(input_dir, output_dir, metadata, reference_db, negative_controls, config, run_id):
    """Run complete ASV generation pipeline."""
    
    pipeline = ASVPipeline(config_file=Path(config) if config else None)
    
    results = pipeline.run_complete_pipeline(
        input_dir=Path(input_dir),
        output_dir=Path(output_dir),
        metadata_file=Path(metadata) if metadata else None,
        reference_db=Path(reference_db) if reference_db else None,
        negative_controls=list(negative_controls) if negative_controls else None,
        run_id=run_id
    )
    
    click.echo("ASV pipeline completed successfully!")
    click.echo(f"Results saved to: {results['output_structure']['run_dir']}")
    click.echo(f"Final ASV FASTA: {results['final_asv_fasta']}")
    click.echo(f"Final ASV table: {results['final_asv_table']}")


if __name__ == '__main__':
    asv_cli()