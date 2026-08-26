# DE RULE 30-LENS · strategie-document (2026-08-26)

Bron: analyse van de Rule 30-rechterrand/linkerhelft (het "toy universe"-
onderzoek: permutiviteit, diagonaal-periodiciteit, alternatieve
geschiedenissen, multiverse-vertakking; incl. Miles Wilsons berekening
tot rij 2^46). Dit document vertaalt die methodologie naar ons
Collatz-programma en definieert kandidaat-taken voor Campagne XVI.
Alle taken zijn patroonherkenning/formalisatie-klasse (geen zware
computes, conform randvoorwaarden 19-aug).

## 1. De Rule 30-feiten, compact

1. **Linker-permutiviteit**: flip het linker-invoerveld in elk van de 8
   gevallen en de uitvoer flipt. Gevolg: informatie in de beginstaat
   blijft langs zijn diagonaal behouden; het tijdomkeer-argument dwingt
   dan dat elke rechterrand-diagonaal periodiek is VANAF rij 0.
2. **Periodes zijn machten van 2**: duiventil op drie kleuren-triples
   plus permutiviteit geeft periode(volgende diagonaal) | 2^(n+1).
3. **Randobservabele**: lengte van het rechtste blauwe blok per rij
   geeft 1,3,1,4,1,3,1,6,1,3,1,4,1,3,1,7,... met de herhalingswet:
   nieuw getal alleen op posities 2^m, en na elk nieuw getal replay van
   de volledige prefix. De innovatiereeks 1,3,4,6,7,9,... is ondanks
   9 biljard bits rekenwerk (rij 2^46) volslagen ondoorzichtig.
4. **Alternatieve geschiedenissen**: achteruit in de tijd is de
   voortzetting alleen uniek ONDER de achtergrond-aanname (zwart blijft
   zwart); een blok zwart kan ook uit blauw komen. Lokale geneste
   structuur: de beginstaat verschijnt oneindig vaak opnieuw in
   zichzelf (elke voldoende grote randdriehoek is een recapitulatie).
5. **Linkerhelft = informatieverlies**: diagonalen uiteindelijk
   periodiek maar niet vanaf de start; vertakkingspunten (waar
   zaad-informatie asymptotisch kan overleven) bestaan alleen direct na
   een effen zwarte diagonaal waarvan de linkerbuur-cyclus een EVEN
   aantal blauwe velden heeft; zeldzaam (eerste bij diagonaal 53208).
   Zoekstrategie: itereer direct in cykel-ruimte, niet in
   beginstaat-ruimte.
6. **Centrale kolom**: drie open vragen (herhaling, gelijkverdeling,
   shortcut) · exact de vorm van onze orbit-vragen.

## 2. Geverifieerd resultaat: het rand-isomorfisme (Obs 610)

Script: `scripts/403_rule30_lens.py` (draait schoon).

De Collatz-klokreeks **B(k) = v2(3^k − 1)** voldoet aan de IDENTIEKE
herhalingswet als Rule 30's randreeks A. Scherpste vorm: beide reeksen
zijn **ruler-functies** · door de herhalingswet hangt de waarde alleen
van v2(k) af:

    A(k) = g(v2(k)),  g: 0,1,2,3,4,5 → 1,3,4,6,7,9   (g ONBEKEND, open)
    B(k) = f(v2(k)),  f(0)=1, f(j)=2+j                (BEWEZEN, LTE)

Ze stemmen exact overeen op de dichte verzameling {v2(k) ≤ 2} (7/8 van
alle posities) en verschillen exact op de diepe kloksloten {v2(k) ≥ 3}.

**Lezing**: Collatz bevat een OPGELOSTE Rule 30. Het 2-adische uurwerk
(Obs-reeks repunit-mechanisme: ord(3 mod 2^j) = 2^(j−2), verdubbelende
periodes) is structureel dezelfde machine als Rule 30's rechterrand,
maar onze innovatiefunctie is affien en bewezen waar die van Rule 30
ondoorgrondelijk is. De hardheid die Rule 30 in zijn randreeks bewaart,
zit bij Collatz verplaatst naar (a) de Sturmian/rotatielaag van log2(3)
en (b) de carry-koppeling tussen de bases · precies onze Conversie-These.

## 3. Het woordenboek

| Rule 30 | Collatz (ons dossier) |
|---|---|
| linker-permutiviteit: info behouden per diagonaal | delings-sweep injectief per digit gegeven carry (Thm 36); backward stairway rule: in-graad staat in de ternaire staart (R216) |
| informatie-verliezende linkerhelft | min/branch-aggregatie, w-keuzes; nul kanaalcapaciteit (Thm 22) |
| diagonaalperiodes = machten van 2 | ord(3 mod 2^j) = 2^(j−2), verdubbelend uurwerk (repunit-mechanisme) |
| randreeks blauwe-bloklengtes | v2(3^k−1) ruler-reeks (Obs 610, bewezen kant) |
| innovatiereeks g (open, gerekend tot 2^46) | innovatie f(j)=2+j (gesloten vorm) · de opgeloste kant |
| alternatieve geschiedenissen (achtergrond kan blauw zijn) | 2-adische realisaties van periodieke streams: −5, −17 (F17/F19); ladder-dekpunt = n=−5 |
| achtergrond-aanname ⟹ unieke geschiedenis | positiviteit/integraliteit ⟹ het vermoeden (realisatie-census: elke samenzwering is een 2-adisch punt, geen enkele positief) |
| lokale geneste structuur: beginstaat recapituleert in zichzelf | families a·3^k−1 = zelfde uurwerk, andere fase; torenwet; monden-structuur |
| lichtkegel: rechtste 10 cellen bepalen 10 diagonalen voorgoed | n mod 2^t bepaalt t pariteiten (Terras); Thm 9 shadowing; klasse-DP i.p.v. integers |
| zoeken in cykel-ruimte i.p.v. beginstaten | shape/census-machinerie (onafhankelijk in beide projecten uitgevonden: methodologische validatie) |
| vertakkingspunten: zeldzaam, alleen na effen zwarte diagonaal met even blauwtelling | needle-survivors / rijke-kaste-runs / spring-blocked rungs (r=+0.52); het monsterprofiel = de configuratie waar geheugenloosheid puntsgewijs zou moeten falen |
| centrale kolom: herhaling? gelijkverdeling? shortcut? | orbit-statistiek: Thm C-mixing (gelijkverdeling); Thm 22 nul-opslag (shortcut); cykel-census (herhaling) |

## 4. Strategie: kandidaat-taken Campagne XVI (Rule 30-lens)

Gerangschikt op kosten; alles binnen patroonherkenning-modus.

**S1 · Permutiviteits-audit van de delings-automaat (formalisatie,
gratis).** Thm 36 geeft de exacte twee-regel-automaat (append-1 +
pariteits-sweep). Stel exact vast in welke richting deze permutief is
(digit ↦ (digit', carry') injectief gegeven carry-in: deling door 2 is
inverteerbaar) en waar de informatie sterft (de w-keuze / min-branch).
Leid vervolgens de uurwerk-periodiciteit met machten-van-2 af als
corollarium in Rule 30-stijl (duiventil + permutiviteit + tijdomkering).
Doel: één verenigde stelling "het 2-adische uurwerk = de permutieve rand
van de Collatz-automaat". Kandidaat voor NOTE-stelling.

**S2 · Familie-innovatiereeksen (goedkoop experiment).** Onze f is
affien voor de repunit-familie (a=1). Families a·3^k−1 zijn zelfde klok
met faseverschuiving a en rotatie-offset log2(a). Vraag: bestaat er een
familie waarvan de innovatiereeks NIET affien is? Zo ja, dan hebben we
"Rule 30-hardheid" exact gelokaliseerd binnen Collatz; zo nee, dan is de
tamheid van de klok universeel over families en zit alle hardheid
aantoonbaar in de rotatielaag. Beide uitkomsten zijn winst.

**S3 · De achtergrond-stelling (formalisatie, gratis).** Rule 30's
"geschiedenis is uniek onder de achtergrond-aanname" is letterlijk onze
situatie: achteruit-continuering via de backward stairway rule is
meerduidig, en de alternatieve geschiedenissen zijn exact de 2-adische
punten (−5, −17, ladder-dekpunten). Schrijf dit als schone sectie:
het vermoeden = "er bestaat geen positieve-integer alternatieve
geschiedenis". Verbindt F17/F19/realisatie-census met één beeld;
kandidaat-sectie voor een paper 7 of de arXiv-note.

**S4 · Vertakkingspunt-jacht op de geordende zijde (goedkoop
experiment).** Rule 30's les: zelfs op de "regelmatige" kant kan
maat-nul zaad-informatie overleven, maar ALLEEN bij zeldzame exacte
configuraties (effen zwart + even telling), en die vind je door direct
in cykel/stream-ruimte te itereren. Onze nul-opslag-stellingen (Thm 20,
22) zijn dichtheids-niveau; het monsterprofiel (Obs 578/595) is de
puntsgewijze uitzonderingskandidaat. Taak: classificeer in stream-ruimte
de exacte configuratie-voorwaarde waaronder continuering NIET geforceerd
is (analoog aan "even blauwtelling"), en census die voorwaarde. Dit is
de Rule 30-vorm van de open kern E★.

**S5 · Tijdomkeer-systematiek (middel).** Rule 30 wint "periodiek vanaf
rij 0" door injectiviteit achterstevoren te draaien. Wij hebben de
achteruit-mint geometrisch(1/3) (R216, fundament dichtheidsrecord).
Taak: inventariseer welke voorwaarts GEMETEN wetten achterstevoren
GEFORCEERD zijn (onvoorwaardelijk), zoals de boom-verdunning. Elke
achteruit-afleiding promoveert een meting naar een stelling.

**S6 · Methodologische adopties (permanent).**
- **Lichtkegel-truncatie**: bereken alleen wat de observabele kan
  beïnvloeden (Wilsons 128-diagonalen-truc = onze klasse-DP; als
  principe vastleggen in RUNBOOK).
- **Shear/coördinaatkeuze**: periodiciteit wordt pas zichtbaar na de
  juiste transformatie (hun shear = onze stairway/macro-coördinaten).
- **Innovatie-epistemologie**: rapporteer expliciet wat berekend vs
  bewezen is; een reeks 9 biljard bits diep zonder patroon is een
  eerlijk en publiceerbaar negatief resultaat (conform onze
  drie-benige standaard).

## 5. Prioriteit

S1 en S3 zijn puur formalisatiewerk op bestaand bewezen materiaal en
kunnen direct; S2 en S4 zijn goedkope experimenten met beslisbare
uitkomst; S5 is een inventarisatie-pass over het dossier. Voorstel
volgorde: S1 → S2 → S3 → S4 → S5.

## 6. Status (bijgewerkt 2026-08-26)

- **S1 KLAAR (Obs 611)**: flip-pariteit oneven op elk niveau j=3..20;
  de klok vertakt nooit · Rule 30's verdubbelingscriterium herbewijst
  LTE in automaat-vorm en de klok-multiverse is bewijsbaar triviaal.
  Script: `scripts/404_family_innovation.py`.
- **S2 KLAAR + BESLIST (Lemma 38, Obs 611)**: familie-klok-dichotomie
  mod 8 (vlak voor a≡5,7; affiene toren voor a≡1,3, bewezen via
  index-2-coset). Geen Rule 30-hardheid in de klok-laag; alle hardheid
  zit bewijsbaar in de koppeling klok×roulette. Bijvangst:
  vlakke-klok-coset = rijke refill-kaste (alle alternators).
- **S3 OPEN**: achtergrond-stelling als papersectie (volgende stap).
- **S4 OPEN**: vertakkingspunt-census in stream-ruimte; na S1 is de
  scherpe vraag: de klok levert géén branch-punten, dus elke
  overlevende zaad-informatie moet via de w-keuze/min-branch lopen ·
  census die configuraties.
- **S5 OPEN**: tijdomkeer-inventarisatie.
