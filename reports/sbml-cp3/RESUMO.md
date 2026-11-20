# Validação SBML · CP3

Rodamos os 1.823 casos semânticos em duas revisões do DNAr. A primeira
reproduziu os 350 aprovados informados pelo Papel B. A segunda incluiu o
detector de cinética do Papel C e manteve os 350 casos, com mais quatro
aprovações.

| Status | CP1 | B | B+C |
| --- | ---: | ---: | ---: |
| Passou | 12 | 350 | 354 |
| Falhou | 269 | 184 | 180 |
| Não suportado | 1.542 | 1.289 | 1.289 |

O CP1 usava `react()` e um filtro de recursos anterior. A comparação que
isola a entrada do C é B → B+C. Nela, as duas rodadas usaram a SBML Test
Suite 3.5.0, R 4.4.1, libSBML 5.21.2 e `react4()` com deSolve/LSODA 1.40
(`rtol=1e-10`, `atol=1e-12`). Foi selecionada uma variante por caso, com
preferência por L3V2. O limite por ponto foi `1e-6 + 1e-4 × abs(referência)`.
Com as tolerâncias dos próprios casos, passaram 349 de 453 trajetórias no B
e 353 de 452 no B+C.

## O que mudou

Passaram os casos **00022, 00804, 00806 e 00807**. No 00804, o maior erro
absoluto caiu de 0,204 para cerca de 2,06e-9. Nenhum dos casos aprovados pelo
B perdeu esse status.

O **01066** já falhava no B, mas ainda gerava trajetória. No B+C, a
importação para no detector de cinética. A comparação
`if (abs(e - g) > 1e-09)` de `sbml_match_stoichiometry()` recebe um valor
ausente ao analisar a estequiometria atribuída. Para reproduzir na revisão
do C, rode `Rscript scripts/run_sbml_baseline.R --cases=01066` com a suíte
baixada.

O **00001** é um exemplo aprovado no B+C. No CP1 o importador parava ao
avaliar o símbolo `compartment` na lei cinética. Com B, a trajetória ficou
a cerca de 1,86e-12 da referência no ponto de maior erro absoluto. O
**00025**, que usa `FunctionDefinition`, continua filtrado antes da simulação.

## Compatibilidade observada

| Recurso no B+C | Passou | Falhou | Não suportado |
| --- | ---: | ---: | ---: |
| `BoundaryCondition` | 72 | 20 | 155 |
| `InitialAssignment` | 86 | 86 | 151 |
| `NonUnityCompartment` | 62 | 19 | 217 |
| `RateRule` | 0 | 0 | 411 |
| `EventNoDelay` | 0 | 0 | 277 |

Um caso pode ter várias tags. A matriz mede o resultado de casos que
contêm cada recurso. A matriz de nível e versão cobre as variantes
escolhidas na execução, com 5 aprovações L2V5, 4 L3V1 e 345 L3V2 no B+C.

## Ambiente e revisão

O `R CMD INSTALL` passou na imagem Docker reconstruída. A correção remove
o erro de sintaxe em `SIGNALS_N_SYSTEMS.R` e deixa os scripts de demonstração
e diagnóstico fora do carregamento do pacote. A cópia local com #1, B e C
combinados também instalou e passou nos 31 testes selecionados de SBML,
testbench, cinética e tradução. A revisão do testbench do B conferiu a
troca para `react4()`, a remoção do filtro de `InitialAssignment` e a
inclusão de saídas de parâmetros e compartimentos.

As rodadas completas medem as propostas de B e C nos commits abaixo. A
rodada final deve ser executada no commit integrado que o grupo congelar.

- B em `032edf2` ([casos](antes/cases.csv), [manifesto](antes/manifest.json)).
- B+C em `58e242b` ([casos](depois/cases.csv), [manifesto](depois/manifest.json)).
- [Transições por status](comparacao/transicoes.csv) e [casos que mudaram](comparacao/mudancas.csv).
- [Matriz por recurso](comparacao/matriz-recursos.csv) e [por nível e versão](comparacao/matriz-nivel-versao.csv).
- [Instruções de execução](../../examples/sbml/README.md).
