<!-- Served as the site's /linking page (rendered in the browser by web/linking.html).
     Relative links below resolve against the site's base path, e.g. `docs` is
     /foundervariants/docs, so they are written for the website, not for a file viewer. -->
# FounderVariants — Linking to the UI

How another site links straight into a FounderVariants search: for example, a Mouse
Phenome Database GWAS page for a Diversity Outbred study, linking from its peak to the
8 founders' genotypes at that locus. The link opens the normal web UI with the search box
and filters filled in and the results already showing. The user can then keep exploring.

Everything here is a **plain GET link**: an `<a href>` with no form, no POST, no
JavaScript on your side and no API key. These parameters are a public contract. Changing
or removing one is a breaking change, recorded in the project changelog.

For the JSON/CSV data API (use it when you want the *data*, not the page), see the
[API reference](docs).

---

## Base URL

```
https://churchilllab.jax.org/foundervariants/
```

Parameters go after the `?`. The trailing `/` is part of the address.
The original address, `https://churchilllab.jax.org/foundersnps/`, still works: it
redirects here with the parameters intact. Use `/foundervariants/` in links you
publish, since it saves the redirect.

## The 30-second version: link to a GWAS peak

```
https://churchilllab.jax.org/foundervariants/?release=GRCm38&chr=10&pos=128350502
```

This opens the GRCm38 database and shows every founder SNP within **50 kb either side** of
chr10:128,350,502. The row at that exact position is highlighted and scrolled into view.
Here that row is rs29358506, a missense SNP.

To show a different window, add `window`:

```
https://churchilllab.jax.org/foundervariants/?release=GRCm38&chr=10&pos=128350502&window=250kb
```

If you build links in code from values your page already has, you need only the assembly,
the chromosome and the peak position:

```js
function founderSnpsLink(assembly, chr, peakBp, windowBp = 50000) {
  const p = new URLSearchParams({release: assembly, chr, pos: peakBp, window: windowBp});
  return "https://churchilllab.jax.org/foundervariants/?" + p;
}
founderSnpsLink("GRCm38", "10", 128350502)
// -> https://churchilllab.jax.org/foundervariants/?release=GRCm38&chr=10&pos=128350502&window=50000
```

### Try it

Each of these opens this site with the link applied:

- [A GWAS peak, default ±50 kb window](?release=GRCm38&chr=10&pos=128350502):
  896 SNPs, the peak (rs29358506, a missense SNP in *Cs*) highlighted
- [The same peak, ±1 Mb, coding changes only](?release=GRCm38&chr=10&pos=128350502&window=1Mb&consequence=missense_variant,stop_gained,start_lost):
  a wide window made readable by a filter
- [A whole gene, with its lead SNP marked](?release=GRCm38&gene=Cs&mark=rs29358506)
- [SNPs private to CAST around the peak](?release=GRCm38&chr=10&pos=128350502&window=250kb&private_to=CAST):
  the peak SNP drops out, and the page says why
- [A GRCm39 range](?release=GRCm39&chr=11&start=100.5Mb&end=100.6Mb&mark=100550004)
- [A link with mistakes in it](?release=GRCh38&gene=Pten&private_to=BALB&consequence=bogus):
  the page still works and lists what it ignored

---

## Parameters

Every parameter is optional. With none, the page opens empty, as it always has.

### `release`: which genome build

| value | opens |
|---|---|
| `1505`, `GRCm38`, `mm10` | Sanger REL-1505 on GRCm38 (Ensembl 75 annotation) |
| `2021`, `GRCm39`, `mm39` | Sanger REL-2112 on GRCm39 (Ensembl 104 annotation) |

Case-insensitive. **Without `release`, the page opens the newest assembly, currently
GRCm39.** So **always pass the assembly your positions are on.** A GRCm38 coordinate
looked up in GRCm39 lands somewhere else, silently: the page shows real SNPs, just not the
ones at your locus. This matters for the MPD DO GWAS files, which are on GRCm38 (the
`.GRCm38.` in their file names): those links need `release=GRCm38`. Assembly names
(`GRCm38`) are the safest choice: they keep working if a newer database for the same
assembly replaces the current one.

### What to search: give one of these

They are checked in this order, and the first one present wins.

| parameter(s) | example | searches |
|---|---|---|
| `q` | `q=19:5.0-5.2Mb` | exactly what a user would type in the search box: a region, gene symbol or rsID |
| `region` | `region=11:100.5-102Mb` | a region |
| `gene` | `gene=Pten` | a gene's full span: a symbol (case-insensitive) or an Ensembl gene ID (`gene=ENSMUSG00000013663`). A symbol naming several loci (e.g. `U6`) shows them to pick from, so link by Ensembl ID when a symbol isn't unique |
| `rsid` | `rsid=rs30847235` | one SNP |
| `chr` + `pos` (+ `window`) | `chr=10&pos=128350502&window=1Mb` | `pos ± window`; **`pos` is also marked** (below). `window` defaults to **50 kb** |
| `chr` + `start` + `end` | `chr=11&start=100.5Mb&end=102Mb` | that range |

- **Positions** are in base pairs. You can also write them with `kb`/`Mb`
  (`128.35Mb`), and commas are allowed (`128,350,502`). The one exception is `mark`,
  described below.
- **`chr`** accepts `10` or `chr10`, and `X`, `Y`, `MT`.
- **`window`** is a half-width: `window=1Mb` searches 1 Mb on each side of `pos`. The window
  stops at position 1 at the start of a chromosome.

### `mark`: highlight specific rows

`mark=128350502` or `mark=rs29358506`, or several separated by commas:
`mark=rs29358506,128400000`. Matching rows are highlighted, and the first one is
scrolled into view. `chr`+`pos` marks `pos` automatically. Use `mark` to highlight
anything else, e.g. the lead SNP when you linked by `gene`.

Because `mark` is comma-separated, **positions inside `mark` must not use thousands
separators**. `mark=1,500,000` means three positions. Write `1500000` or `1.5Mb`.

If no row matches, the page says why. Either the position isn't a founder SNP or a filter
removed it, or it's past the row limit (see [Row limit and window size](#row-limit-and-window-size)).

### Filters

These set the matching UI controls, which the user can then change.

| parameter | example | effect |
|---|---|---|
| `consequence` | `consequence=missense_variant,stop_gained` | only variants with *any* of these consequences. Names are the Sequence Ontology terms listed by [`GET /{release}/so_terms`](2021/so_terms) |
| `private_to` | `private_to=CAST` | only SNPs where this founder differs from the other seven. Accepts short or full names: `A/J`, `B6`, `129`, `NOD`, `NZO`, `CAST`, `PWK`, `WSB`, or `CAST_EiJ` etc. |
| `sdpn` | `sdpn=43658` | one exact strain-distribution pattern: a number encoding which founders carry which allele, shown in the results' SDPN column |
| `low_confidence` | `low_confidence=true` | also include sites where a founder's genotype call is low-confidence (`true`/`1`/`yes`) |
| `limit` | `limit=10000` | maximum rows; default **2000**, maximum **25000** |

`differ` (two founders that carry different bases) is available in the [API](docs) but
not in the UI, because the page has no control for it.

---

## Row limit and window size

The page shows at most `limit` rows (default 2,000), **in position order**. Founder SNPs
are dense. Measured around the example peak on GRCm38:

| window around chr10:128,350,502 | high-confidence founder SNPs |
|---|---|
| ±50 kb (the default) | 896 |
| ±100 kb | 2,068 |
| ±250 kb | 5,499 |
| ±1 Mb | 21,111 |

So a ±1 Mb window with the default limit shows only the first 2,000 rows, i.e. the
**left edge** of the window, and the peak is cut off. The page detects this and says so,
but it's better to avoid:
- **Keep the window tight** (the 50 kb default), or
- **raise `limit`** (up to 25,000 rows; the page stays usable at that size), or
- **filter**, e.g. `consequence=missense_variant,stop_gained` turns tens of thousands of
  rows into a handful.

Density varies along the genome, roughly with how divergent the founders are locally, so
these numbers show the scale but don't predict every locus.

## What the page does with a link

- **Bad parameters never break the page.** An unknown release, founder or consequence, or
  a malformed position, is ignored. A one-line note under the search box says which
  parameter and why. The rest of the link still applies.
- **The address bar follows the user.** After the link loads, every search or filter
  change rewrites the URL in a short, standard form (e.g. `chr`+`pos` becomes `q=` plus
  `mark=`). Copying the address bar always reproduces what's on screen, so users can share
  any search, whether or not it started from a link. This uses `history.replaceState`, so
  live typing doesn't fill the Back button history.
- **Marks belong to the link.** Once the user types a different search, the highlights and
  any link notes clear.

## Why GET and not POST

A link-out sends a build, a locus and a few filters, under 200 characters. Browsers and
the server accept URLs many times that long. A GET link can be bookmarked, emailed, opened
in a new tab and refreshed without a "resubmit form?" prompt. POST would only be worth it
for sending a large list, e.g. hundreds of SNPs or loci at once. The page doesn't support
that: it shows one search at a time, so a list would need a batch feature, not a different
way of arriving.

## When you want the data, not the page

Link to the API instead. Its selectors are `region`, `gene` and `rsid`, without `q` or
`chr`/`pos`, and it can return a file:

```
https://churchilllab.jax.org/foundervariants/1505/variants?region=10:128300502-128400502&format=csv&limit=5000
```

That downloads one row per SNP with a column for each founder. See the [API reference](docs).
