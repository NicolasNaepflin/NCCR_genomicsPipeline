"""Annotation of external proteins through the existing CLI entry point."""
import re
import shutil

import pytest
import yaml
from click.testing import CliRunner

from workflow.main import main


@pytest.mark.integration
@pytest.mark.parametrize('external', [True, False])
def test_annotate_protein_inputs(repo_root, tmp_path, capfd, external):
    if shutil.which('snakemake') is None:
        pytest.skip('Snakemake not installed')
    out = tmp_path / 'output'
    config = {'outDir': str(out), 'eggnog_db': str(tmp_path / 'database')}
    if external:
        # No read samplesheet; spaces exercise shell quoting.
        proteins = tmp_path / 'existing proteins.faa'
        proteins.write_text('>gene1\nMKKLLV\n')
        config['preexisting_genes'] = {'genome_a': str(proteins),
                                       'genome_b': str(proteins)}
    else:
        proteins = out / 'assembly' / 'genome_a' / 'genome_a.faa'
        proteins.parent.mkdir(parents=True)
        proteins.write_text('>gene1\nMKKLLV\n')
        config['sample'] = ['genome_a']
    config_file = tmp_path / 'config.yaml'
    config_file.write_text(yaml.safe_dump(config))

    result = CliRunner().invoke(main, ['annotate', '-c', str(config_file), '--dry'])
    captured = capfd.readouterr()
    output = result.output + captured.out + captured.err
    assert result.exit_code == 0, output
    expected = 'emapper_genes' if external else 'emapper'
    assert re.search(rf'^rule {expected}:', output, re.MULTILINE), output
    assert not re.search(r'^rule (assemble_wga|prokka|qc):', output, re.MULTILINE)
    if external:
        assert 'genome_a.emapper.annotations' in output
        assert 'genome_b.emapper.annotations' in output
        assert f"'{proteins}'" in output
        assert re.search(r'^total\s+3$', output, re.MULTILINE), output


@pytest.mark.integration
@pytest.mark.parametrize('genes', [{}, [], {'bad/name': 'genes.faa'}])
def test_invalid_preexisting_genes(tmp_path, capfd, genes):
    if shutil.which('snakemake') is None:
        pytest.skip('Snakemake not installed')
    config_file = tmp_path / 'config.yaml'
    config_file.write_text(yaml.safe_dump({
        'outDir': str(tmp_path / 'out'), 'eggnog_db': '/db',
        'preexisting_genes': genes,
    }))
    result = CliRunner().invoke(main, ['annotate', '-c', str(config_file), '--dry'])
    captured = capfd.readouterr()
    assert result.exit_code != 0
    assert 'preexisting_genes' in result.output + captured.out + captured.err
