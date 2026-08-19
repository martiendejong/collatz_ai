# arXiv-pakket — status en checklist

**STATUS: KLAAR VOOR VERZENDING · WACHT OP EXPLICIETE GO VAN MARTIEN.**
Er wordt NIETS verzonden zonder zijn expliciete opdracht (vaste regel,
zie NOTE.md Obs 545-addendum en de Lagarias-afspraak van 25 juli 2026).

## Inhoud van het pakket

1. **density_one.tex / .pdf** (33 pagina's, compileert schoon, draft 0.2)
   - Hoofdclaim: assemblage van het gamma -> 1-programma met alle
     exponenten bewezen; vijf boekhoudtaken expliciet; certificaten t/m
     k = 21 (gamma = 0.9184, 3.49 miljard constraints, exact geheeltallig).
   - Augustus-secties: vezeltoren + eindige inductie (dominantie op het
     eindpunt exact-rationaal gecertificeerd, k = 13/14); Doeblin-lemma
     (Q < 1 gecertificeerd met worst-case staartopslag); klok-camouflage-
     stellingen; Booleaans-spectraal instrument (het regenererende bos).
   - Statuslabels [PROVED]/[MEASURED]/[BOOKKEEPING] op elke bewering;
     vier eerlijke retractaties gedocumenteerd in het publieke log.

2. **Aanbevolen categorieen**: math.NT (primair), math.DS (secundair).

3. **Voor verzending nog te doen (mechanisch, ~1 uur, op go):**
   - [ ] Bibliografie aanvullen (Tao 2019, Terras 1976, Barina 2025,
         Hercher 2023, Halbeisen-Hungerbuhler 1997 - deels al aanwezig)
   - [ ] Log-referenties bevriezen op een vaste commit-hash
   - [ ] Laatste spellingscontrole beide talen van citaten
   - [ ] arXiv-account/endorsement regelen (math.NT vergt mogelijk
         endorsement voor een nieuw account)

## Companion-materiaal (niet voor arXiv, wel openbaar)
- De cursus: https://martiendejong.github.io/collatz_ai/ (28 lessen)
- Het volledige onderzoekslog: research/NOTE.md (Obs 1-600+)
- Alle scripts (397+) en certificaten in de repo

## Waarom publiceren de juiste zet is (advies, 19 aug 2026)
Het programma-knelpunt is verschoven van ontdekking naar bewijstechnologie
die wij niet bezitten. De reductie is af: een open lemma, alles eromheen
stelling/certificaat. De Booleaans-spectrale sectie spreekt een vakgebied
aan (analyse van Booleaanse functies) waarvan het gereedschap nog nooit op
Collatz is gericht. Publicatie dateert het werk en nodigt de juiste
specialisten uit op onze scherpste formulering.
