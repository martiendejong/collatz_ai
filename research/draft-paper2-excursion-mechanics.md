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

## Novelty-check TODO (vóór indiening)
- [ ] Renewal-wet u·2^s − v: de v=1-vorm (2^k a −1 → 3^k a −1) is folklore;
      is de algemene vorm mét 3x−1-subtrahend-interpretatie ergens gepubliceerd?
      (zoek: "Collatz" + "3x−1 dual", Böhm–Sontacchi-parametrisaties)
- [ ] Mod-8-selectieregels / LTE-dood-stop: zoek in Wirsching, Lagarias-
      annotated-bibliography (2 delen), Applegate–Lagarias inverse tree.
- [ ] n²-recordwet-afleiding: Kontorovich–Lagarias (2009) "stochastic models"
      geven Brownian/LD-heuristieken voor stopping times — hebben zij de
      excursie-exponent c=1 / m*=4/3 expliciet? (vermoedelijk κ-niveau wel,
      de exacte identiteit + compositiewet vermoedelijk niet)
- [ ] Schaduw-decompositie langs negatieve cycli: geen referentie bekend;
      meest waarschijnlijk-nieuwe onderdeel samen met de camouflage-stelling.
- [ ] Aperiodiciteits-dichotomie: p-adische repulsie is klassiek; is deze
      toepassing (eventueel-periodieke itinerary ⟹ integer-cyclus) al ergens
      genoteerd? (zoek: Bernstein–Lagarias 2-adische conjugatie)

## Positionering
- Paper 1 (dichtheidsrecord x^0.902, draft-arxiv-note.md) = het harde
  onvoorwaardelijke resultaat; dit paper = de structuurlaag. Kunnen apart;
  paper 1 eerst (record verjaart).
- Alles elementair verifieerbaar (geen zware machinerie): geschikt voor
  bijv. Experimental Mathematics / J. Number Theory / Acta Arithmetica
  (advies Lagarias vragen in dezelfde mail).
- Volledige reproduceerbaarheid: scripts 177–181 in de publieke repo.
