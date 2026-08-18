import re, json

# ---- 1. b7l4: begrensde-wiebel-stelling in subject 3 ----
s = open('b7l4.html', encoding='utf-8').read()
anchor = 'Everything except\n      "and that cannot continue forever" is a theorem.</span>\n    </div>'
assert anchor in s, "anker b7l4 subject 3 niet gevonden"
add = anchor + '''
    <p class="nl">En sinds augustus 2026 is dit keurslijf nog een maat strakker: de
    <b>begrensde-wiebel-stelling</b>. Elk aperiodiek patroonwoord met begrensde discrepantie
    (elk Sturmian-woord, van welke helling ook) is voor positieve gehele getallen onrealiseerbaar ·
    bewijs in drie regels: begrensde discrepantie betekent een begrensde baan, een begrensde baan
    van een geheel getal bezoekt eindig veel waarden en wordt dus periodiek, en een periodieke baan
    draagt een periodiek woord · tegenspraak. Het monster moet dus niet alleen aperiodiek en
    recurrent zijn, maar ook <b>onbegrensd van zijn eigen gemiddelde wegdrijven</b>: op elke schaal
    steeds langere excursies in zijn eigen patroon. Elke nette wiskundige constructie is daarmee
    uitgesloten; alleen het zelf-versnellende woord rest.</p>
    <p class="en">And since August 2026 this straitjacket is one size tighter: the
    <b>bounded-wobble theorem</b>. Every aperiodic pattern word of bounded discrepancy (every
    Sturmian word, of any slope) is unrealisable by positive integers · a three-line proof: bounded
    discrepancy means a bounded orbit, a bounded integer orbit visits finitely many values and thus
    becomes periodic, and a periodic orbit carries a periodic word · contradiction. So the monster
    must be not only aperiodic and recurrent but must also <b>drift unboundedly from its own
    average</b>: ever longer excursions within its own pattern at every scale. Every tidy
    mathematical construction is thereby excluded; only the self-accelerating word remains.</p>'''
s = s.replace(anchor, add, 1)
open('b7l4.html', 'w', encoding='utf-8').write(s)

# ---- 2. b7l2: vijf lagen -> zes lagen ----
s = open('b7l2.html', encoding='utf-8').read()
o1 = 'camouflage van de afbeelding op vijf lagen een stelling: symbolen, kanalen, telling, kern én\n      klok.'
n1 = ('camouflage van de afbeelding op zes lagen een stelling of exacte meting: symbolen, kanalen, telling, kern, '
      'klok · en de koppeling zelf: de gezamenlijke verdeling van kop-hoek (Weyl-rotatie) en staart-klok is gemeten '
      'volmaakt onafhankelijk (χ² = 228 bij 225 vrijheidsgraden over 480.000 stappen).')
assert o1 in s
s = s.replace(o1, n1, 1)
o2 = "map's camouflage is a theorem on five layers: symbols, channels, counting, kernel and\n      clock."
n2 = ("map's camouflage is a theorem or exact measurement on six layers: symbols, channels, counting, kernel, "
      'clock · and the coupling itself: the joint distribution of head angle (Weyl rotation) and tail clock is '
      'measured perfectly independent (χ² = 228 at 225 degrees of freedom over 480,000 steps).')
assert o2 in s
s = s.replace(o2, n2, 1)
open('b7l2.html', 'w', encoding='utf-8').write(s)

# ---- 3. b5l5: k=18-datapunten in het drie-instrumenten-subject ----
s = open('b5l5.html', encoding='utf-8').read()
o3 = 'creep dies (DENSITY) or persists (CEILING) is Conjecture G.'
n3 = ('creep dies (DENSITY) or persists (CEILING) is Conjecture G. Latest depth (k = 18): the TR instrument reads '
      '0.8357 and the γ̄-rate instrument 0.8404 · both still rising, same band, the fork still open.')
assert o3 in s
s = s.replace(o3, n3, 1)
o4 = 'is Conjecture G genoemd'
if o4 not in s:
    # NL variant: zoek de NL-zin over de vork/kruip in subject 4
    o4 = 'DENSITY) of blijft (CEILING) is Conjecture G.'
    n4 = ('DENSITY) of blijft (CEILING) is Conjecture G. Laatste diepte (k = 18): het TR-instrument meet 0.8357 en '
          'het γ̄-rate-instrument 0.8404 · beide nog stijgend, zelfde band, de vork nog open.')
    if o4 in s:
        s = s.replace(o4, n4, 1)
open('b5l5.html', 'w', encoding='utf-8').write(s)

# ---- 4. b2l4: 100M rotaties vermelden ----
s = open('b2l4.html', encoding='utf-8').read()
o5 = 'patroon plus alle 22 rotaties zijn al in seconden weerlegd.'
n5 = ('patroon plus alle 22 rotaties zijn al in seconden weerlegd; van ticket 1 zijn inmiddels de eerste '
      '100 miljoen rotaties afgelopen (O(1) per rotatie via de baan van het rationale vaste punt) · nul treffers.')
assert o5 in s
s = s.replace(o5, n5, 1)
o6 = 'pattern plus all 22 rotations are already refuted in seconds.'
n6 = ('pattern plus all 22 rotations are already refuted in seconds; for ticket 1 the first 100 million rotations '
      'have now been swept (O(1) per rotation via the orbit of the rational fixed point) · zero hits.')
assert o6 in s
s = s.replace(o6, n6, 1)
open('b2l4.html', 'w', encoding='utf-8').write(s)

# validatie
for fn in ['b7l4.html', 'b7l2.html', 'b5l5.html', 'b2l4.html']:
    t = open(fn, encoding='utf-8').read()
    o = len(re.findall(r'<div[ >]', t)); c = t.count('</div>')
    assert o == c, (fn, o, c)
    for blob in re.findall(r'<script type="application/ld\\+json">(.*?)</script>', t, re.S):
        json.loads(blob)
print('vier lessen bijgewerkt en gevalideerd')
