#!/bin/bash

# Verifica se o nome do arquivo foi passado como argumento
if [[ $# -eq 0 ]]; then
    echo "Uso: $0 <nome_do_arquivo>"
    exit 1
fi

# Lê o nome do arquivo fornecido como argumento
arquivo=$1

# Verifica se o arquivo existe
if [[ ! -f $arquivo ]]; then
    echo "O arquivo $arquivo não existe."
    exit 1
fi

# Cria um arquivo temporário para armazenar o código sem comentários
temp=$(mktemp)

# Remove os comentários do arquivo, preservando a identação do código
sed 's/\s*!.*//' "$arquivo" > "$temp"

# Substitui o conteúdo do arquivo original pelo conteúdo do arquivo temporário
mv "$temp" "$arquivo"

# Informa que a remoção de comentários foi concluída com sucesso
echo "Comentários removidos com sucesso do arquivo $arquivo."

