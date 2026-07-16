# Het Collatz-onderzoek in 200 woorden

Neem een getal. Is het even: halveer het. Is het oneven: keer drie plus een.
Herhaal. Het vermoeden (Collatz, 1937): je komt altijd op 1 uit. Niemand
heeft het ooit kunnen bewijzen.

Wij onderzochten het met een AI, in ~900 onderzoeksronden. Wat we vonden:

- Een wereldrecord: bewezen dat minstens x^0,91 van de getallen onder x
  op 1 uitkomt (het oude record stond sinds 2003 op 0,84).
- De machine erachter is een uurwerk: vermenigvuldigen met 3 is volledig
  voorspelbaar; al het raadsel zit in het delen door 2.
- Het proces gedraagt zich als een perfect eerlijke munt: geen geheugen,
  geen patroon dat een getal kan uitbuiten - exact geteld, niet geschat.
- Waarom een piano nooit zuiver stemt en waarom Collatz moeilijk is, is
  letterlijk hetzelfde feit: machten van 2 en 3 passen nooit precies op
  elkaar.
- De afdaling naar 1 wordt aangedreven door de bovenkant van het getal -
  iets wat alleen echte, eindige getallen hebben.

Opgelost is het niet: bewijzen dat elk getal daalt vergt wiskunde die nog
niet bestaat. Maar de kaart is scherper dan ooit: zes artikelen, alles
openbaar op github.com/martiendejong/collatz_ai.
