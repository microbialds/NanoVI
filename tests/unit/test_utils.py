"""Unit tests for general utility functions (bin/utils.py).

Covers gather_fastq_paths (filesystem traversal) and timed_function (decorator).
Functions that shell out to external tools (fastplong, minimap2) are excluded
and covered by nf-test integration tests.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'bin'))
from utils import gather_fastq_paths, timed_function


class TestGatherFastqPaths:
    def test_finds_fastq_files_in_directory(self, tmp_path):
        (tmp_path / "sample.fastq").write_text("@r\nACGT\n+\nIIII")
        result = gather_fastq_paths([str(tmp_path)])
        assert any("sample.fastq" in p for p in result)

    def test_finds_fq_extension(self, tmp_path):
        (tmp_path / "sample.fq").write_text("")
        result = gather_fastq_paths([str(tmp_path)])
        assert any(".fq" in p for p in result)

    def test_finds_fastq_gz_files(self, tmp_path):
        (tmp_path / "sample.fastq.gz").write_bytes(b"")
        result = gather_fastq_paths([str(tmp_path)])
        assert any(".fastq.gz" in p for p in result)

    def test_finds_fq_gz_files(self, tmp_path):
        (tmp_path / "sample.fq.gz").write_bytes(b"")
        result = gather_fastq_paths([str(tmp_path)])
        assert any(".fq.gz" in p for p in result)

    def test_ignores_non_fastq_files(self, tmp_path):
        (tmp_path / "notes.txt").write_text("not a fastq")
        (tmp_path / "data.csv").write_text("a,b")
        result = gather_fastq_paths([str(tmp_path)])
        assert result == []

    def test_returns_sorted_paths(self, tmp_path):
        for name in ["c.fastq", "a.fastq", "b.fastq"]:
            (tmp_path / name).write_text("")
        result = gather_fastq_paths([str(tmp_path)])
        assert result == sorted(result)

    def test_recurses_into_subdirectories(self, tmp_path):
        sub = tmp_path / "subdir"
        sub.mkdir()
        (sub / "nested.fastq").write_text("")
        result = gather_fastq_paths([str(tmp_path)])
        assert any("nested.fastq" in p for p in result)

    def test_direct_file_path_accepted(self, tmp_path):
        fq = tmp_path / "direct.fq"
        fq.write_text("")
        result = gather_fastq_paths([str(fq)])
        assert any("direct.fq" in p for p in result)

    def test_empty_input_returns_empty_list(self):
        assert gather_fastq_paths([]) == []

    def test_returns_absolute_paths(self, tmp_path):
        (tmp_path / "sample.fastq").write_text("")
        result = gather_fastq_paths([str(tmp_path)])
        assert all(os.path.isabs(p) for p in result)

    def test_multiple_directories_combined(self, tmp_path):
        dir_a = tmp_path / "a"
        dir_b = tmp_path / "b"
        dir_a.mkdir()
        dir_b.mkdir()
        (dir_a / "s1.fastq").write_text("")
        (dir_b / "s2.fastq").write_text("")
        result = gather_fastq_paths([str(dir_a), str(dir_b)])
        assert len(result) == 2


class TestTimedFunction:
    def test_returns_function_result(self):
        result = timed_function("double", lambda x: x * 2, 21)
        assert result == 42

    def test_passes_positional_args(self):
        result = timed_function("add", lambda a, b: a + b, 3, 4)
        assert result == 7

    def test_passes_keyword_args(self):
        def add(a, b=0):
            return a + b
        result = timed_function("add_kw", add, 3, b=7)
        assert result == 10

    def test_propagates_exception(self):
        with pytest.raises(ZeroDivisionError):
            timed_function("bad", lambda: 1 / 0)
