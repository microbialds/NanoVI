import sys
import types
import importlib
import math
from contextlib import contextmanager

class FakeArray(list):
    def __truediv__(self, other):
        return [x / other for x in self]

def setup_stubs():
    fake_numpy = types.ModuleType("numpy")
    fake_numpy.array = lambda x: FakeArray(x)
    fake_numpy.exp = lambda x: [math.exp(v) for v in x]
    sys.modules['numpy'] = fake_numpy

    fake_flat = types.ModuleType("flatten_dict")
    fake_flat.unflatten = lambda d: d
    sys.modules['flatten_dict'] = fake_flat

    bio_seq = types.ModuleType("Bio.Seq")
    bio_seq.Seq = lambda s: s
    sys.modules['Bio.Seq'] = bio_seq

class FakeAlignment:
    def __init__(self, query_name, ref_name, stats):
        self.query_name = query_name
        self.reference_name = ref_name
        self._stats = stats
    def get_cigar_stats(self):
        return [self._stats]

def make_pysam(alignments):
    fake_pysam = types.ModuleType('pysam')
    class FakeAlignmentFile:
        def __init__(self, path, mode='rb'):
            self.alignments = alignments
        def fetch(self):
            return self.alignments
    fake_pysam.AlignmentFile = FakeAlignmentFile
    return fake_pysam

class FakePool:
    def __init__(self, threads):
        self.threads = threads
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
    def map(self, func, chunks):
        return list(map(func, chunks))

@contextmanager
def patched_env(alignments):
    setup_stubs()
    sys.modules['pysam'] = make_pysam(alignments)
    import multiprocessing
    original_pool = multiprocessing.Pool
    multiprocessing.Pool = FakePool
    try:
        yield
    finally:
        multiprocessing.Pool = original_pool

import pathlib

def load_vi():
    vi_path = pathlib.Path(__file__).parents[1] / 'bin' / 'variational_inference.py'
    spec = importlib.util.spec_from_file_location('variational_inference', vi_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_zero_alignments():
    with patched_env([]):
        vi = load_vi()
        log_p, zeros, longest = vi.get_cigar_op_log_probabilities('dummy', threads=4)
        assert log_p == []
        assert zeros == [0, 1, 2, 3]
        assert longest == {}

def test_few_alignments():
    alignments = [
        FakeAlignment('q1', '1:ref', [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]),
        FakeAlignment('q2', '1:ref', [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    ]
    with patched_env(alignments):
        vi = load_vi()
        log_p, zeros, longest = vi.get_cigar_op_log_probabilities('dummy', threads=4)
        assert len(log_p) == len(vi.CIGAR_OPS) - len(zeros)
        assert 'q1' in longest and 'q2' in longest

def test_threads_zero_equivalent_to_one():
    alignments = [
        FakeAlignment('q1', '1:ref', [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    ]
    with patched_env(alignments):
        vi = load_vi()
        result_zero = vi.get_cigar_op_log_probabilities('dummy', threads=0)
        result_one = vi.get_cigar_op_log_probabilities('dummy', threads=1)
        assert result_zero == result_one
