# Fluxo de Potência para Sistemas de n Barras

Este repositório implementa um código para o cálculo do fluxo de potência em sistemas elétricos de potência com n barras. O usuário deve definir os parâmetros do sistema no arquivo `dados.py`.

## Estrutura do Repositório

- `.gitignore`: Arquivo que especifica os arquivos e diretórios a serem ignorados pelo Git.
- `Trabalho_Sistemas_de_Potencia.pdf`: Documento com informações teóricas ou explicações sobre o trabalho.
- `dados.py`: Arquivo onde o usuário deve definir os dados do sistema elétrico (número de barras, parâmetros das linhas, cargas, etc.).
- `exemplo1.m`: Arquivo de exemplo em MATLAB.
- `gauss-seidel.py`: Implementação do método de Gauss-Seidel para resolver o fluxo de potência.
- `newton-raphson-desacoplado.py`: Implementação do método de Newton-Raphson desacoplado.
- `newton-raphson.py`: Implementação do método de Newton-Raphson para resolver o fluxo de potência.

## Requisitos

- Python 3.8 ou superior
- Bibliotecas Python necessárias:
  - `numpy`
  - `scipy`

Para instalar as dependências, execute:
```bash
pip install -r requirements.txt
```

## Como Utilizar

1. Edite o arquivo `dados.py` para configurar o sistema elétrico com base no número de barras, dados das linhas, geração e carga.
2. Escolha o método desejado para resolver o fluxo de potência:
   - Para o método de Gauss-Seidel, execute:
     ```bash
     python gauss-seidel.py
     ```
   - Para o método de Newton-Raphson, execute:
     ```bash
     python newton-raphson.py
     ```
   - Para o método de Newton-Raphson desacoplado, execute:
     ```bash
     python newton-raphson-desacoplado.py
     ```
3. Os resultados do fluxo de potência serão exibidos no terminal.


