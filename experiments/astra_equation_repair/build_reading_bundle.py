"""Build a single public Markdown reading file from the scientific artifacts.

No fitting, evaluation, or alteration of scientific results occurs here.
"""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
REPOSITORY = 'joelpaulson/ChiralML'
BRANCH = 'astra-equation-repair'
SUBDIR = 'experiments/astra_equation_repair'
BASE = f'https://github.com/{REPOSITORY}/blob/{BRANCH}/{SUBDIR}/'
RAW = f'https://raw.githubusercontent.com/{REPOSITORY}/{BRANCH}/{SUBDIR}/'


def link(target, image=False):
    if target.startswith(('http://', 'https://', '#')):
        return target
    return (RAW if image else BASE) + target


def include(path):
    text = path.read_text()
    text = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', lambda m: f'![{m[1]}]({link(m[2], True)})', text)
    text = re.sub(r'(?<!!)\[([^\]]*)\]\(([^)]+)\)', lambda m: f'[{m[1]}]({link(m[2])})', text)
    return text


def main():
    sections = ['REPORT.md', 'FINAL_EQUATION.md', 'EXPERIMENT_LOG.md']
    manifest = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in sections}
    intro = '''# ChiralML equation-repair experiment: complete reading document

This document combines the scientific report, full fitted equation, chronological
proposal/evaluation log, and numerical verification for one exploratory study.
It is generated from the original artifacts without changing the scientific results.

**Conclusion:** retain the original five-term B32 structure. Three reasoning-driven
edits failed to provide convincing transfer improvement. The final coefficients
use all 773 observations; the final chemistry-holdout assessment used only the
earlier development-fitted coefficients. Discovery leakage and imperfect blindness
to one published literature example are explicitly documented below.

All following sections are study content. The complete code, workbook, predictions,
figures and source provenance are available through the linked repository artifacts.

'''
    parts = [intro, f'[Experiment directory](https://github.com/{REPOSITORY}/tree/{BRANCH}/{SUBDIR}) · '
             f'[Reproduction instructions]({BASE}README.md) · [Raw reading document]({RAW}CHATGPT_READING.md)\n\n']
    for name in sections:
        parts += ['\n---\n\n', f'Source: [{name}]({BASE}{name})\n\n', include(HERE/name), '\n']
    parts += ['\n---\n\n# Numerical verification\n\n```json\n', (HERE/'verification.json').read_text(), '```\n\n']
    parts += ['# Artifact index\n\nThe source workbook and symbolic-regression implementation are also present on this branch.\n\n',
              f'- [Source workbook](https://github.com/{REPOSITORY}/blob/{BRANCH}/data_analysis/data/analysis_data.xlsx)\n',
              f'- [Existing SyMANTIC implementation](https://github.com/{REPOSITORY}/tree/{BRANCH}/symbolic_regression)\n']
    for p in sorted(HERE.rglob('*')):
        rel = p.relative_to(HERE)
        if not p.is_file() or any(part.startswith('.') or part == '__pycache__' for part in rel.parts):
            continue
        if p.name in ('CHATGPT_READING.md', 'reading_sources.json'):
            continue
        parts.append(f'- [{rel.as_posix()}]({BASE}{rel.as_posix()})\n')
    (HERE/'CHATGPT_READING.md').write_text(''.join(parts))
    (HERE/'reading_sources.json').write_text(json.dumps(dict(source_sha256=manifest, repository=REPOSITORY, branch=BRANCH), indent=2)+'\n')
    print(f'Wrote CHATGPT_READING.md ({(HERE/"CHATGPT_READING.md").stat().st_size:,} bytes)')


if __name__ == '__main__':
    main()
