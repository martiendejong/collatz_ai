# Paper 2 (skelet): The Arithmetic Mechanics of Collatz Excursions

**Status: SKELET voor publicatiereview 2026-07-25. Alle resultaten bewezen +
machinegeverifieerd (scripts 177–181); novelty-check per onderdeel nodig vóór claims.**

## Voorgestelde titel
*The arithmetic mechanics of Collatz excursions: binomial normal form,
negative-cycle shadowing, and the n² record law*

## Kernboodschap (abstract-concept)
Elke klim van een Collatz-baan is exacte algebra, geen toeval: we geven de
gesloten klim-motor T(u·2^s − v) = 3u·2^{s−e} − (3v−1)/2^e (de subtrahend
volgt de 3x−1-map, d.w.z. het 2-adische schaduwdoel doorloopt zijn eigen
negatieve baan), klassificeren alle klimmechanismen (LTE-exitwetten,
Mersenne-dood-stop, mod-8-selectieregels), decoderen de beroemde baan van 27
knoop voor knoop, en leiden uit één exacte identiteit —
I(4/3) = log₂3 − 4/3 met I′(m) = log₂(2(m−1)/m) —
de decennia-oude empirische n²-wet van padrecords af, inclusief de
compositievoorspelling P(e=1) = 3/4 voor lange klimmen (gemeten: 0.714/0.725).
Complementair bewijzen we dat de symboollaag exact i.i.d. Geom(1/2)² is en
re-entry-dieptes exact geometrisch: alle arithmetiek zit in de residu-draad,
en klim-"geluk" heeft geen statistische structuur.

## Secties

1. **Setup & notation** — T, macro-stap n = 2^K a − 1, itinerary coding,
   coding-injectiviteit (lemma, éénregel).
2. **The binomial normal form** — renewal-stelling (tweeregel-bewijs),
   Mersenne = vast punt v=1, eindige-brandstof, 703-golf als voorbeeld
   (2^15 − 3^5 → 3·2^12 − 91 → 9·2^8 − 17).
3. **Exit laws and selection rules** — K′ = v₂(3^K a + 1) − 1 (l=1);
   LTE op de Mersenne-ruggengraat; dood-stop-stelling; kanaaltabel mod 8;
   ±-uitsluiting als motor van de klim-correlatie.
4. **The orbit of 27, decoded** — de volledige keten met regelverificatie
   op elke stap (en de vier andere <10^6-records idem, machinegecheckt).
5. **Negative-cycle shadowing** — schaduwstelling (depth-verlies exact e per
   stap), cyclusvolgorde-tracking (meetdata), per-episode-tekort 0.415d,
   schaduwdiepte als lokale klimvoorspeller (data: kanteling bij d ≈ 7).
6. **Statistical camouflage** — exacte onafhankelijkheid (l, K′) ~ Geom²;
   re-entry exact geometrisch (9 octaven). Consequentie: κ is onverbeterbaar
   op symboolniveau.
7. **The n² law** — de identiteit, het minimum bij m* = 4/3, kosten
   1 bit/bit, P(exc ≥ E) ≍ 2^{−E}, piek ~ n², de 3/4-compositiewet +
   metingen; aansluiting bij Roosendaal-records en Kontorovich–Lagarias-
   stochastische modellen (onze bijdrage: de afleiding + compositie).
8. **The aperiodicity dichotomy** — stelling + bewijs; gevolg: periodieke
   ruggengraat van de divergentieverzameling bevat geen integers.
9. **Outlook** — complexiteitsladder (automatisch/Mahler als volgende trap),
   3-adische geschiedeniscodering, relatie tot Tao 2019 (aggregaat vs per-baan).

## Novelty-check — UITGEVOERD 2026-07-24 (Obs 347)
- [x] **n²-recordwet: GEDEGRADEERD.** Lagarias–Weiss 1992 (Ann. Appl. Prob. 2,
      229–261) voorspelden lim log t(n)/log n = 2 al via RRW-model + large
      deviations. §7 hergeformuleerd: "herafleiding binnen het episode-
      raamwerk"; ónze bijdragen: de gesloten identiteit I(4/3) = log₂3 − 4/3
      met I′(m) = log₂(2(m−1)/m), de compositiewet Geom(3/4), en de eerste
      empirische compositieverificatie op records (0.7775/0.1854).
      RESTCHECK: L–W integraal lezen — staat de 4/3-tilt er expliciet in?
- [x] **Aperiodiciteits-dichotomie: BEKEND** (Terras/Everett-codering;
      Bernstein–Lagarias 1996, Canad. J. Math. 48, 1154–1169). §8 wordt
      expositie mét attributie; ons 2-adische-repulsie-bewijs als alternatief
      bewijs presenteren, niet als nieuwe stelling.
- [x] **Normaalvorm/schaduw/selectieregels/camouflage: GEEN treffers** —
      dit is de kandidaat-nieuwe kern. v=1-klim (2^k−1 → 3^k−1) is folklore:
      als zodanig attribueren. Te citeren ankers: Kontorovich–Sinai 2003
      (2-adische equidistributie, bij §6), Kontorovich–Lagarias
      arXiv:0910.1944, 3x+d-cyclusliteratuur (Cox arXiv:2101.04067,
      arXiv:2101.08060) bij de negatieve-cyclusfamilie.
- [ ] RESTCHECK vóór indiening: Wirsching (LNM 1681) + Lagarias' annotated
      bibliography (2 delen) scannen op (u,s,v)-vorm en selectieregels;
      Lagarias–Weiss 1992 integraal.

## Bonus voor Paper 1 (uit dezelfde check)
Chunlei Liu, "Counting the Collatz numbers" (arXiv:2512.13760, dec 2025)
bewijst x^0.3227 en noemt expliciet "The historical record is 0.84" —
onafhankelijke bevestiging (dec 2025!) dat ons x^0.902 een record is.
Citeren in draft-arxiv-note.md.

## Positionering
- Paper 1 (dichtheidsrecord x^0.902, draft-arxiv-note.md) = het harde
  onvoorwaardelijke resultaat; dit paper = de structuurlaag. Kunnen apart;
  paper 1 eerst (record verjaart).
- Alles elementair verifieerbaar (geen zware machinerie): geschikt voor
  bijv. Experimental Mathematics / J. Number Theory / Acta Arithmetica
  (advies Lagarias vragen in dezelfde mail).
- Volledige reproduceerbaarheid: scripts 177–181 in de publieke repo.
