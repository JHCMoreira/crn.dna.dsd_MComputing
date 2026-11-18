# Resultado inicial dos testes SBML

Rodada sobre o importador da base `acdd21a`, com a SBML Test Suite 3.5.0.
Foram contabilizados 1.823 casos semânticos, sem IDs repetidos. Foi usada
uma variante por caso, com preferência por L3V2.

| Resultado | Casos |
| --- | ---: |
| Passou | 12 |
| Falhou | 269 |
| Não suportado | 1.542 |

Os limites usados foram `1e-6` absoluto e `1e-4` relativo, combinados por
ponto. Com as tolerâncias originais dos casos, passaram 7 das 48 trajetórias
que chegaram à comparação. Esses dois critérios não devem ser confundidos.

Os casos não suportados foram excluídos pelo filtro de recursos ausentes no
importador. Os demais podem falhar durante a importação, na simulação ou na
comparação. A coluna `reason` da [tabela completa](cases.csv) registra o motivo.

## Casos analisados

**00048 — aprovado.** Reação S1 → S2, com k=1, S1 inicial de 0,015 e S2
inicial zero, em um compartimento de dimensão zero. As referências são
S1(t)=0,015 exp(-t) e S2(t)=0,015(1-exp(-t)). Nos 51 instantes entre 0 e 5,
o maior erro absoluto foi 1,17e-12. Os CSVs de
[simulação](example-00048-simulation.csv) e
[referência](example-00048-reference.csv) permitem conferir a comparação.

**00025 — não suportado.** Usa `FunctionDefinition`, que o importador atual
não extrai nem expande na fórmula cinética. Foi excluído pelo filtro,
sem simulação ou validação nessa rota.

**00001 — falha de importação.** A lei cinética contém
`compartment * k1 * S1`. O tamanho do compartimento é lido para as condições
iniciais, mas seu símbolo não entra no ambiente de avaliação da fórmula.
O erro é `object 'compartment' not found`. Esse caso básico é uma prioridade
para a correção da leitura de compartimentos pelo Papel B.

## Leitura das tabelas

- [Por componente](component_tags.csv) e [por tag](test_tags.csv): as contagens
  se sobrepõem, pois um caso pode ter várias tags.
- [Causas](causes.csv): triagem por tags e mensagem; exige revisão manual.
- [Nível/versão](compatibility.csv): descreve apenas as variantes selecionadas.
- [Ambiente](manifest.json) e [versões do R](sessionInfo.txt): parâmetros da rodada.

Validação vazia significa não avaliada. Quando preenchida, corresponde à
função `validate_sbml()` existente, que registra erros de leitura do libSBML
sem chamar `checkConsistency()`. Uma aprovação vale para aquele caso e
tolerância; não certifica suporte completo ao SBML.

A próxima rodada deve usar as mesmas entradas e parâmetros após os merges
de B/C, revendo o filtro de recursos que passaram a ser suportados.

## Referências

- [Método e execução](../../examples/sbml/README.md)
- Casos oficiais [00048](https://github.com/sbmlteam/sbml-test-suite/blob/3.5.0/cases/semantic/00048/00048-model.m),
  [00025](https://github.com/sbmlteam/sbml-test-suite/blob/3.5.0/cases/semantic/00025/00025-model.m) e
  [00001](https://github.com/sbmlteam/sbml-test-suite/blob/3.5.0/cases/semantic/00001/00001-sbml-l3v2.xml).
