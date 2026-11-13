#!/usr/bin/env python3
"""Assemble the isolated v2 draft; original manuscript is never modified.

This is the initial assembly. Subsequent evidence updates are made to paper/v2/main.tex.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'paper/v2'
if (DEST / 'main.tex').exists():
    raise SystemExit('The v2 manuscript already exists. This historical assembly '
                     'script must not overwrite subsequent scientific revisions.')
text = (ROOT / 'paper/main.tex').read_text()

def section(start, end, filename):
    global text
    a, b = text.index(start), text.index(end)
    text = text[:a] + (DEST / filename).read_text().strip() + '\n\n' + text[b:]

def substitute(old, new):
    global text
    if text.count(old) != 1:
        raise ValueError(f'Expected unique replacement anchor: {old[:100]!r}')
    text = text.replace(old, new, 1)

a = text.index(r'\begin{abstract}')
b = text.index(r'\end{abstract}', a) + len(r'\end{abstract}')
text = text[:a] + (DEST / 'abstract_revision.tex').read_text().strip() + text[b:]
section(r'\section{Introduction}', r'\section{Model}', 'intro_revision.tex')
section(r'\section{Theory near the voter point}', r'\section{Discussion}', 'theory_revision.tex')
section(r'\section{Discussion}', r'\begin{acknowledgments}', 'discussion_revision.tex')

substitute('every rule with $w(j)>0$ for', 'every rule with $w(0)=0$ and $w(j)>0$ for')
substitute('Moreover, since also $w(j)<1$, every configuration of a finite lattice reaches consensus with',
           'For $0<p<1$ and $0<\\alpha<1$ in the parameter range studied, $w(j)<1$ as well. Every configuration of a finite lattice then reaches consensus with')
substitute('They have no basin of attraction: a\nsingle defect destroys them, and random initial conditions never end in them for $q>0$. Still,',
           'They were not observed as final states of the random-start simulations at $q>0$, and prepared patterns were typically destroyed by a single defect (Appendix~\\ref{app:frozen}). Still,')
substitute('which covers all the $q$ studied\nhere.', 'which covers the principal $\\alpha=1/2$ family studied here; it does not cover the $\\alpha=1$ checks.')
substitute('at a DP or Ising-type transition it vanishes continuously.',
           'at a DP absorption transition it vanishes continuously. An Ising transition between two active phases instead concerns the magnetization and need not make $\\rho$ vanish.')
substitute('There is no window in which\nboth droplets shrink: the two phases exchange stability at a single point.',
           'No interval in which both droplets shrink is resolved by these measurements. Together with the orientation comparison, this supports a common coexistence point within the measured resolution.')
substitute('25 of 32 runs at $L = 128$ decayed', '25 of 32 runs at $L = 128$ reached an absorbing state')
substitute('none of the 16 runs at $L = 256$ or of the 8 runs each at $L = 512$\nand 1024 decayed within',
           'none of the 16 runs at $L = 256$ or of the 8 runs each at $L = 512$\nand 1024 reached an absorbing state within')
substitute('Below $\\qs$ the lifetime\nsaturates in $L$; above $\\qs$ it grows faster than any power of $L$ and leaves the accessible\nwindow between $L = 128$ and 256, as expected for a phase that can only be destroyed by a\nfluctuation of the size of the system.',
           'The absorption counts differ from crossings of the half-density threshold: at $L=128$ the latter counts are 30 of 32 and 11 of 32 at $1.008\\,\\qs$ and $1.018\\,\\qs$, respectively. Below $\\qs$ the half-density lifetime saturates in $L$; above it the lifetime increases rapidly and exceeds the observation window at larger sizes. These censored observations support the stability of the active phase but do not determine an asymptotic growth law in $L$.')
substitute('In a DP- or GV-type transition, by contrast, a single seed survives with finite probability\neverywhere in the active phase, with no size threshold.',
           'Here $R_c$ is the radius at which the finite-time survival probability is one half, not a strict minimum size: smaller seeds have a nonzero chance to grow. The useful nucleation diagnostic is the growth of this probabilistic radius and its relation to the independently measured front velocity.')
substitute('At a GV\npoint the same holds, and at a DP point $\\delta_{\\rm eff}\\to 0.451$.',
           'GV criticality has voter-like asymptotic spreading, but this does not impose a finite-time bound on $\\delta_{\\rm eff}$ in a nonlinear model~\\cite{castellano2012,azizi2020}; at a DP point $\\delta_{\\rm eff}\\to0.451$.')
substitute('At $(0.99, 0.01)$ the decay is still accelerating at\n$t = 10^6$ ($\\delta_{\\rm eff} = 2.1$, $P_s(10^6) = 1.0\\times10^{-6}$ against $1.9\\times10^{-5}$ on the\nvoter curve),',
           'At $(0.99,0.01)$, $\\delta_{\\rm eff}$ reaches about 2.1 near $t=2.5\\times10^5$. At $t=10^6$, only three of $3\\times10^6$ seeds remain, giving $P_s=10^{-6}$ with a 95\\% binomial interval $[0.21,2.92]\\times10^{-6}$; the late-time tail does not determine a precise exponent,')
substitute('Together, these results place a first-order transition at every $p\\le 0.995$ we studied: a GV\npoint would need a $q_c$ at which the bulk is active, seeds decay no faster than voter seeds, and\nwalls broaden diffusively, and none of the points between $\\qs$ and $1.2\\,\\qs$ at\n$p = 0.98$--$0.995$ behaves like that.',
           'These complementary measurements support continuation of the first-order line toward the voter point. Their interpretation is least secure at $p=0.995$, where threshold uncertainty and crossover times are largest. Accelerated seed extinction is supporting evidence for a nucleation barrier, rather than a universal finite-time test excluding GV criticality. We distinguish this finite-resolution continuation from the well-resolved coexistence at $p=1/2$ and from any claim about every fixed $p<1$.')
substitute('Near the voter point the threshold is exactly\nproportional to $\\alpha$:', 'Near the voter point the measured threshold is consistent with proportionality to $\\alpha$:')
substitute('The first-order\ntransition is not an artifact of synchronous updating.',
           'These tests show that the coexistence signatures survive substantial dilution of synchronous updates. The strict continuous-time limit requires an extrapolation in $r$; the measured finite-$r$ results already exclude a dependence on fully synchronous updating alone.')
substitute('In the voter limit the closure conserves $m$ and gives $\\rho^* = (z-2)(1-m^2)/[2(z-1)]$.',
           'In the voter limit the continuous-time closure conserves $m$ and gives $\\rho^*=(z-2)(1-m^2)/[2(z-1)]$. The synchronous map differs: at $m=0$ its active fixed point is $\\rho^*=1/2$.')
substitute('The active phase has no metastable branch.',
           'No metastable active branch is resolved by these runs.')
substitute('For this rule the transition is thus first order on one side only: the absorbing state is\nmetastable where it loses, the active phase is not.',
           'The two protocols therefore reveal a marked asymmetry: small seeds are suppressed where a macroscopic active region invades, while no stationary retreating active slab is resolved. Since the density data do not determine the limiting jump, we leave the order of this transition open.')
substitute('The band of absorbing-state metastability,\n$1.5^\\circ$--$2^\\circ$ wide, is Case~4A of Ref.~\\cite{cox2013} in two dimensions. It is what the\nsearch of Ref.~\\cite{molofsky1999} missed: it is narrower than the spacing of their parameter\nscans, and inside it their seeds of one or two sites lie below $R_c$.',
           'The band of suppressed initial seed growth is about $1.5^\\circ$--$2^\\circ$ wide. It resembles the invasion-versus-small-seed contrast of Case~4A in Ref.~\\cite{cox2013}, but finite-time seed decay does not establish asymptotic linear stability. Its narrowness provides a possible explanation for why the earlier scan did not resolve such behavior.')

text = '% Version 2: independent analysis and additional campaigns, November 2026.\n' + text
(DEST / 'main.tex').write_text(text)
(DEST / 'refs.bib').write_text((ROOT / 'paper/refs.bib').read_text() + '\n' + (DEST / 'refs_additions.bib').read_text())
labels = set(re.findall(r'\\label\{([^}]+)\}', text))
refs = set(re.findall(r'\\(?:eq)?ref\{([^}]+)\}', text))
print('Wrote', DEST / 'main.tex')
print('Undefined internal references:', sorted(refs-labels))
