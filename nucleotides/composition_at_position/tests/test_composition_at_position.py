from click.testing import CliRunner
import filecmp

from nucleotideCountPerPosition import cli


def test1_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--column', 'A',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 1


def test2_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--column', 'A', ])
    assert result.exit_code == 2


def test3_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--column', 1,
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest1.tsv")


def test4_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--column', 2,
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest2.tsv")


def test5_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input2.tsv',
                                 '--column', "A",
                                 '--header',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest1.tsv")


def test6_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input2.tsv',
                                 '--column', "B",
                                 '--header',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest2.tsv")

def test7_nucleotideCountPerPosition():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input2.tsv',
                                 '--column', "B",
                                 '--header',
                                 '--chunksize', 1,
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest2.tsv")
