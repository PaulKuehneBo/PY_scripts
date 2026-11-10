Ich zeichne die Regeln direkt aus dem aktuellen Kraftmodell und Integrator nach: welche Paarungen Beiträge erzeugen, welche davon auf Atomkerne übertragen werden und welche Kräfte tatsächlich in die Bewegung eingehen. Danach fasse ich getrennt zusammen, was sich bei der vorgeschlagenen Beschränkung auf Elektronen freier Paare ändern würde.


Die Berechnung ist teilweise bereits partikelbasiert: `CoulombForceModel.calculate` erzeugt einzelne Kraftbeiträge zwischen `ParticleState`s. Der Integrator fasst diese Beiträge anschließend aber für die Bewegung der Kerne und Elektronen zusammen (`_integrate_electrons`).

## Aktuelle Regeln

| Betroffenes Partikel | Kraftquelle | Was das Modell macht |
|---|---|---|
| **Nucleus** | Andere Nucleonen | Jedes Kernpaar wechselwirkt; beide Kerne erhalten entgegengesetzte Beiträge. Der Integrator verwendet dafür `NUCLEUS_REPULSION_FACTOR`. |
| **Nucleus** | Elektron desselben Moleküls, **gebunden** | Keine direkte Kern-Elektron-Kraft. Wirken auf das gebundene Elektron andere Kräfte, wird dessen gesamte Resultierende als Reaktionskraft auf seinen eigenen Kern übertragen. |
| **Nucleus** | **Ungebundenes** Elektron desselben Moleküls | Auf einen *anderen* Kern wirkt eine Kraft; der Beitrag wird derzeit **nur am Kern** registriert, nicht entgegengesetzt am Elektron. Der eigene Kern des Elektrons ist ausgenommen. |
| **Nucleus und Elektron** | Elektron und Kern aus verschiedenen Molekülen | Beide erhalten entgegengesetzte Beiträge – auch wenn das Elektron zu einer Bindung gehört. Die Kraft auf ein gebundenes Elektron wird zusätzlich auf dessen eigenen Kern übertragen. |
| **Elektron** | Andere Elektronen | Alle Elektronenpaare wechselwirken mit entgegengesetzten Beiträgen, auch die zwei Elektronen einer Bindung oder eines freien Paars. |
| **Elektron** | Partner in einer Bindung oder einem freien Elektronenpaar | Zusätzlich wirkt eine Paarlängen-Constraint. Sie gilt für beide Paararten, nicht für einzelne freie Elektronen. |

Damit ist auch erklärt, warum du bei einem Elektron oft keine Kern-Kraftpfeile siehst: Bei einem gebundenen Elektron im selben Molekül wird die direkte Elektron-Kern-Wechselwirkung übersprungen. Nucleus-Kräfte erscheinen für Elektronen insbesondere bei Wechselwirkungen mit **anderen Molekülen**.

Eine wichtige Asymmetrie: Bei einem ungebundenen Elektron und einem anderen Kern **im selben Molekül** erhält aktuell nur der Kern einen Kraftbeitrag. Das ist keine wechselseitige Paarwechselwirkung. Ein Lauf des Kraftmodells bestätigte genau diese Fallunterscheidung.

## Was „nur freie Elektronenpaare dürfen Kerne beeinflussen“ bedeuten würde

Dafür sollte „freies Elektronenpaar“ eindeutig über `pair_id` **vorhanden** und `bond_id` **nicht vorhanden** erkannt werden. Die Partikelzustände erhalten ihre Paar-IDs in `_particle_states`.

Die Änderung hätte drei Folgen:

1. Direkte Kern-Elektron-Kräfte gäbe es nur noch für Elektronen aus Lone Pairs; einzelne freie Elektronen und Bindungselektronen würden keine Kerne mehr direkt beeinflussen.
2. Wenn auch die bisherige Reaktionskraftübertragung gebundener Elektronen wegfallen soll, würden Kräfte auf Bindungselektronen nicht mehr auf ihren Kern übertragen. Ihre Bewegung könnte den Kern dann nicht mehr über diese Reaktionskräfte mitbewegen.
3. Für Konsistenz sollte die Lone-Pair-Wechselwirkung **auf Kern und Elektron wechselseitig** angewendet werden. Eine Kraft nur auf den Kern würde die bereits beobachtete Asymmetrie fortsetzen.

Dabei muss zusätzlich festgelegt werden, ob Kräfte, die **andere Elektronen** auf ein gebundenes Elektron ausüben und heute über die Reaktionskraft beim Kern ankommen, weiterhin übertragen werden dürfen. Das ist nicht dasselbe wie eine direkte Kernkraft durch ein Bindungselektron.

Noch eine wichtige Unterscheidung beim Debuggen: Die angezeigten Pfeile beruhen auf einem Kraftbericht am aktuellen Visualisierungszustand. Die Bewegung entsteht dagegen aus mehreren Kraftauswertungen während der Integrationsschritte, zusätzlich zu Dämpfung und Brownscher Bewegung. Deshalb müssen die Pfeile nicht exakt die gesamte Kraftgeschichte eines Frames darstellen.

Ich habe hierfür noch keinen Code geändert; das ist die Bestandsaufnahme und die nötige Designentscheidung vor der Umstellung.

## Aktualisierte Regeln und Laufzeit-Schalter

Die Kraftberechnung registriert die Coulomb-Kraft eines ungebundenen Elektrons und eines fremden Nucleus jetzt symmetrisch: Beide Partikel erhalten entgegengesetzte Beiträge. Die Wechselwirkung mit dem eigenen Nucleus bleibt ausgeschlossen.

Kräfte zwischen gebundenen Elektronen und Nucleonen anderer Moleküle sind standardmäßig aktiviert. Mit `b` lassen sie sich während der Simulation umschalten. Die Einstellung gilt für Kraftberichte und alle Kraftauswertungen des Integrators. Ist sie ausgeschaltet, entfallen sowohl die Kraft auf das gebundene Elektron als auch die entgegengesetzte Kraft auf den fremden Nucleus. Die Reaktionskraft, die die verbleibenden Kräfte eines gebundenen Elektrons auf seinen eigenen Nucleus überträgt, bleibt davon unberührt.

Diese Regel gilt nur für die direkte intermolekulare Elektron-Nucleus-Wechselwirkung. Andere Elektron-Elektron-Kräfte sowie Nucleus-Nucleus-Kräfte bleiben aktiv.

## Kubische Abschwächung intermolekularer Elektron-Nucleus-Kräfte

Bei intermolekularen Elektron-Nucleus-Wechselwirkungen basiert die Abschwächung auf dem Coulomb-Abstand `delta` zwischen dem Elektron und dem Ziel-Nucleus. Der Cutoff-Radius ist:

`Bond.minimum_length + Zielatom.electron_radius`

Für ein Elektron, das zu einer Bindung gehört, wird deren `minimum_length` verwendet. Bei einem Lone Pair oder freien Einzelelektron ohne `Bond` gilt der Standardwert von `Bond.minimum_length`. Der Radius stammt vom Zielatom, zum Beispiel `Proton.electron_radius`. Bei `delta >= Cutoff` bleibt die Coulomb-Kraft unverändert; bei `delta < Cutoff` wird sie mit `(delta / Cutoff) ** 3` skaliert und ist bei `delta == 0` gleich null. Da dieser Faktor vor der Erzeugung der entgegengesetzten Beiträge angewendet wird, bleibt die Elektron-Nucleus-Wechselwirkung symmetrisch.

Die bisherige lineare Radius-Abschwächung bleibt für nicht-intermolekulare Elektron-Nucleus-Wechselwirkungen unverändert. Der Laufzeitschalter `b` für gebundene Elektronen wirkt weiterhin zusätzlich und unabhängig von dieser Abschwächung.