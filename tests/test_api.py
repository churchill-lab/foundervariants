"""The HTTP surface, including the discovery behaviour deployment depends on."""
from __future__ import annotations

import csv
import io
import json

import pytest

from foundervariants import encode as enc


class TestService:
    def test_health_lists_the_releases_found(self, client):
        r = client.get("/health")
        assert r.status_code == 200
        assert r.json()["status"] == "ok"
        assert sorted(r.json()["releases"]) == ["8888", "9999"]

    def test_health_touches_no_database(self, client, db_dir):
        # /health is the container healthcheck, so it must not depend on a readable DB.
        import os
        for f in os.listdir(db_dir):
            os.chmod(os.path.join(db_dir, f), 0o000)
        try:
            assert client.get("/health").status_code == 200
        finally:
            for f in os.listdir(db_dir):
                os.chmod(os.path.join(db_dir, f), 0o644)

    def test_genomes_returns_provenance_per_release(self, client):
        data = client.get("/genomes").json()["data"]
        assert len(data) == 2
        assert {d["release"] for d in data} == {"8888", "9999"}
        assert all(d["grc"] and d["source_vcf"] and d["created"] for d in data)

    def test_meta_reports_a_live_count(self, client):
        d = client.get("/9999/meta").json()["data"]
        assert d["n_variants"] == 7
        assert d["founders"] == enc.FOUNDERS

    def test_so_terms(self, client):
        terms = client.get("/9999/so_terms").json()["data"]
        assert [t["term"] for t in terms][0] == "stop_gained"   # severity-ordered
        assert all("bit" in t for t in terms)

    def test_unknown_release_is_404_and_names_what_exists(self, client):
        r = client.get("/1234/meta")
        assert r.status_code == 404
        assert "1234" in r.json()["detail"]

    def test_releases_are_independent(self, client):
        for rel in ("8888", "9999"):
            assert client.get(f"/{rel}/variants?rsid=rs001").json()["meta"]["release"] == rel


class TestGenesEndpoint:
    def test_autocomplete(self, client):
        d = client.get("/9999/genes?q=Pt").json()["data"]
        assert [g["symbol"] for g in d] == ["Pten", "Ptgs2"]

    def test_limit(self, client):
        assert len(client.get("/9999/genes?q=Pt&limit=1").json()["data"]) == 1

    def test_q_is_required(self, client):
        assert client.get("/9999/genes").status_code == 422


class TestVariantsEndpoint:
    def test_region(self, client):
        d = client.get("/9999/variants?region=19:1000-3000").json()["data"]
        assert [v["variants_id"] for v in d["variants"]] == ["rs001", "rs002"]
        assert d["count"] == 2

    def test_region_with_unit_suffix_scales_both_ends(self, client):
        # Regression guard for the parse_region fix, through the HTTP layer.
        a = client.get("/9999/variants?region=19:0.001-0.003Mb").json()["data"]
        b = client.get("/9999/variants?region=19:1000-3000").json()["data"]
        assert [v["variants_id"] for v in a["variants"]] == \
               [v["variants_id"] for v in b["variants"]]

    def test_gene(self, client):
        d = client.get("/9999/variants?gene=Pten").json()["data"]
        assert [v["variants_id"] for v in d["variants"]] == ["rs001", "rs002"]

    def test_rsid(self, client):
        d = client.get("/9999/variants?rsid=rs005").json()["data"]
        assert d["variants"][0]["founders"]["WSB_EiJ"] == "A"

    def test_consequence_multi_term(self, client):
        d = client.get("/9999/variants?region=19:0-99999"
                       "&consequence=stop_gained,missense_variant").json()["data"]
        assert {v["variants_id"] for v in d["variants"]} == {"rs001", "rs002", "rs005"}

    def test_low_confidence_flag(self, client):
        without = client.get("/9999/variants?gene=Pten").json()["data"]["count"]
        with_ = client.get("/9999/variants?gene=Pten&low_confidence=true").json()["data"]["count"]
        assert with_ == without + 1

    def test_founder_patterns(self, client):
        d = client.get("/9999/variants?region=19:0-99999&private_to=WSB_EiJ").json()["data"]
        assert [v["variants_id"] for v in d["variants"]] == ["rs005"]
        d = client.get("/9999/variants?rsid=rs004&differ=A_J,NOD_ShiLtJ").json()["data"]
        assert d["count"] == 1

    def test_exact_sdpn(self, client):
        target = enc.bases_to_sdpn(list("TTTTTTTA"))
        d = client.get(f"/9999/variants?region=19:0-99999&sdpn={target}").json()["data"]
        assert [v["variants_id"] for v in d["variants"]] == ["rs005"]

    def test_no_selector_is_422(self, client):
        assert client.get("/9999/variants").status_code == 422

    def test_bad_region_is_422(self, client):
        assert client.get("/9999/variants?region=garbage").status_code == 422

    def test_unknown_gene_is_404(self, client):
        assert client.get("/9999/variants?gene=NoSuchGene").status_code == 404

    def test_unknown_founder_is_422(self, client):
        r = client.get("/9999/variants?region=19:0-99999&private_to=Nope_J")
        assert r.status_code == 422

    def test_unknown_consequence_is_422(self, client):
        r = client.get("/9999/variants?region=19:0-99999&consequence=nope")
        assert r.status_code == 422

    def test_limit(self, client):
        d = client.get("/9999/variants?region=19:0-99999&limit=1").json()["data"]
        assert d["count"] == 1


class TestExport:
    def _rows(self, body, delim):
        return list(csv.DictReader(io.StringIO(body), delimiter=delim))

    def test_csv_is_flat_with_one_column_per_founder(self, client):
        r = client.get("/9999/variants?gene=Pten&format=csv")
        assert r.status_code == 200
        assert r.headers["content-type"].startswith("text/csv")
        assert "attachment" in r.headers["content-disposition"]
        rows = self._rows(r.text, ",")
        assert len(rows) == 2
        for f in enc.FOUNDERS:
            assert f"founder_{f}" in rows[0]
        assert rows[0]["founder_PWK_PhJ"] == "A"

    def test_tsv(self, client):
        r = client.get("/9999/variants?gene=Pten&format=tsv")
        assert r.headers["content-type"].startswith("text/tab-separated-values")
        assert len(self._rows(r.text, "\t")) == 2

    def test_json_format_is_a_bare_array(self, client):
        r = client.get("/9999/variants?gene=Pten&format=json")
        body = json.loads(r.text)
        assert isinstance(body, list) and len(body) == 2

    def test_unknown_format_falls_back_to_the_envelope(self, client):
        body = client.get("/9999/variants?gene=Pten&format=xlsx").json()
        assert "data" in body and "meta" in body

    def test_empty_result_still_exports(self, client):
        r = client.get("/9999/variants?region=19:500000-600000&format=csv")
        assert r.status_code == 200


class TestStaticAndDownload:
    def test_ui_is_served(self, client):
        r = client.get("/")
        assert r.status_code == 200 and "text/html" in r.headers["content-type"]
        assert 'href="/docs"' in r.text and "%%BASE%%" not in r.text   # no base path: root

    def test_deep_link_serves_the_same_page(self, client):
        # the query string is read client-side (see /linking); the server must not
        # reject or vary on it
        r = client.get("/?release=GRCm38&chr=10&pos=128,350,502&window=1Mb&mark=rs1&junk=%3Cx%3E")
        assert r.status_code == 200 and "parseLink" in r.text
        assert r.text == client.get("/").text

    def test_linking_page_renders_the_vendored_guide(self, client):
        page = client.get("/linking")
        assert page.status_code == 200 and "text/html" in page.headers["content-type"]
        assert "/static/marked.min.js" in page.text and "/static/LINKING.md" in page.text
        assert "://" not in page.text.split("marked.min.js")[0].split("src=")[-1]  # no CDN
        md = client.get("/static/LINKING.md")
        assert md.status_code == 200 and md.text.startswith("<!--")
        assert "# FounderVariants — Linking to the UI" in md.text
        js = client.get("/static/marked.min.js")
        assert js.status_code == 200 and "marked v12.0.2" in js.text[:200]

    def test_main_page_links_to_the_guide(self, client):
        assert 'href="/linking"' in client.get("/").text

    def test_sdp_page_and_its_assets(self, client):
        page = client.get("/sdp")
        assert page.status_code == 200 and "text/html" in page.headers["content-type"]
        for asset in ("/static/site.css", "/static/marked.min.js", "/static/sdpcalc.js"):
            assert asset in page.text
            assert client.get(asset).status_code == 200, asset
        md = client.get("/static/SDP.md")
        assert md.status_code == 200 and "## How SDP and SDPN are calculated" in md.text
        main = client.get("/").text
        assert 'href="/sdp"' in main and "/sdp?sdpn=" in main   # top bar + per-row links
        # the same calculator component in the search page's pop-up and on /sdp
        for asset in ("/static/sdpcalc.js", "/static/sdpcalc-ui.js", "/static/sdpcalc.css"):
            assert asset in main and asset in page.text, asset
            assert client.get(asset).status_code == 200, asset
        assert 'id="sdpdlg"' in main and 'id="sdpcalcbtn"' in main

    def test_linking_page_uses_the_shared_stylesheet(self, client):
        assert 'href="/static/site.css"' in client.get("/linking").text

    def test_guide_has_no_repo_relative_doc_links(self, client):
        # served as web pages: links like (API.md) would 404 there
        import re
        for doc in ("/static/LINKING.md", "/static/SDP.md"):
            assert not re.findall(r"\]\([A-Za-z_-]+\.md", client.get(doc).text), doc

    def test_docs_page_references_the_local_scalar_bundle(self, client):
        body = client.get("/docs").text
        assert "/static/scalar.js" in body        # vendored, never a CDN
        assert "://" not in body.split("scalar.js")[0].split("src=")[-1]

    def test_openapi_schema_is_generated(self, client):
        paths = client.get("/openapi.json").json()["paths"]
        assert "/{release}/variants" in paths

    def test_rqtl2_download_streams_a_sqlite_file(self, client):
        r = client.get("/9999/download/rqtl2")
        assert r.status_code == 200
        assert r.headers["content-type"] == "application/x-sqlite3"
        assert "fv.9999.snps.sqlite" in r.headers["content-disposition"]
        assert r.content[:16] == b"SQLite format 3\x00"


class TestBehindAPathStrippingProxy:
    """FV_BASE_PATH set, as on the server: Apache strips the prefix, so the app receives
    bare paths (/static/x) while every URL it writes carries the prefix."""

    @pytest.fixture
    def proxied(self, db_dir, monkeypatch):
        from fastapi.testclient import TestClient

        from foundervariants import api
        monkeypatch.setenv("FV_DB_DIR", db_dir)
        monkeypatch.setenv("FV_BASE_PATH", "/foundervariants")
        with TestClient(api.create_app()) as c:
            yield c

    def test_static_assets_are_served_at_the_stripped_path(self, proxied):
        # Regression: a Starlette Mount only matched /foundervariants/static/..., so the
        # API docs page's scalar.js 404'd behind Apache
        r = proxied.get("/static/scalar.js")
        assert r.status_code == 200 and "javascript" in r.headers["content-type"]

    def test_pages_carry_the_prefix(self, proxied):
        page = proxied.get("/").text
        assert 'const BASE = "/foundervariants"' in page
        # absolute, so it works from /foundervariants AND /foundervariants/ (a relative
        # href="docs" resolved to /docs from the slash-less URL)
        assert 'href="/foundervariants/docs"' in page
        assert 'href="/foundervariants/linking"' in page
        assert "%%BASE%%" not in page
        guide = proxied.get("/linking").text
        assert 'const BASE = "/foundervariants"' in guide
        assert "src=\"/foundervariants/static/marked.min.js\"" in guide
        assert "%%BASE%%" not in guide
        assert proxied.get("/static/LINKING.md").status_code == 200
        calc = proxied.get("/sdp").text
        assert 'const BASE = "/foundervariants"' in calc and "%%BASE%%" not in calc
        assert 'src="/foundervariants/static/sdpcalc.js"' in calc
        assert 'href="/foundervariants/sdp"' in page
        docs = proxied.get("/docs").text
        assert "src='/foundervariants/static/scalar.js'" in docs
        assert "data-url='/foundervariants/openapi.json'" in docs

    def test_api_routes_answer_at_the_stripped_path(self, proxied):
        assert proxied.get("/health").status_code == 200
        assert proxied.get("/9999/variants?region=1:1-10Mb").status_code == 200

    @pytest.mark.parametrize("bad", ["..%2Fapi.py", "%2E%2E%2Fapi.py", "nope.js", "%2E%2E"])
    # (a bare "../x" is normalised away by the HTTP client before sending; the
    # percent-encoded forms are what actually reach the route)
    def test_static_never_escapes_its_directory(self, proxied, bad):
        assert proxied.get(f"/static/{bad}").status_code == 404


class TestDiscovery:
    def test_a_corrupt_database_is_skipped_not_fatal(self, tmp_path, monkeypatch, db_factory):
        """Deployment relies on this: a half-copied file must degrade, not down, the service."""
        from fastapi.testclient import TestClient

        from foundervariants import api

        d = tmp_path / "dbs"
        d.mkdir()
        db_factory(str(d / "fv.9999.snps.sqlite"))
        (d / "fv.7777.snps.sqlite").write_bytes(b"not a database at all")

        monkeypatch.setenv("FV_DB_DIR", str(d))
        with TestClient(api.create_app()) as c:
            assert c.get("/health").json()["releases"] == ["9999"]

    def test_an_empty_directory_serves_zero_releases(self, tmp_path, monkeypatch):
        from fastapi.testclient import TestClient

        from foundervariants import api
        (tmp_path / "empty").mkdir()
        monkeypatch.setenv("FV_DB_DIR", str(tmp_path / "empty"))
        with TestClient(api.create_app()) as c:
            assert c.get("/health").json()["releases"] == []
            assert c.get("/9999/meta").status_code == 404

    def test_only_the_expected_glob_is_picked_up(self, tmp_path, monkeypatch, db_factory):
        from fastapi.testclient import TestClient

        from foundervariants import api

        d = tmp_path / "dbs"
        d.mkdir()
        db_factory(str(d / "fv.9999.snps.sqlite"))
        db_factory(str(d / "fv.1111.indels.sqlite"))   # wrong infix
        db_factory(str(d / "other.2222.snps.sqlite"))  # wrong prefix
        monkeypatch.setenv("FV_DB_DIR", str(d))
        with TestClient(api.create_app()) as c:
            assert c.get("/health").json()["releases"] == ["9999"]
