$inputFile = "combined_sequences.fasta"
$outputFile = "random_50_sequences_per_cluster.fasta"

# Leer todas las secuencias del archivo de entrada
$sequences = Get-Content -Path $inputFile

# Agrupar secuencias por cluster
$clusters = @{}
$currentCluster = ""
$currentSequence = ""

foreach ($line in $sequences) {
    if ($line -like ">*") {
        if ($currentSequence -ne "") {
            if ($clusters[$currentCluster] -eq $null) {
                $clusters[$currentCluster] = @()
            }
            $clusters[$currentCluster] += ,$currentSequence
            $currentSequence = ""
        }
        $currentCluster = $line.Split('_')[0].TrimStart('>')
        $currentSequence = $line + "`n"
    } else {
        $currentSequence += $line + "`n"
    }
}

# Añadir la última secuencia al cluster correspondiente
if ($currentSequence -ne "") {
    if ($clusters[$currentCluster] -eq $null) {
        $clusters[$currentCluster] = @()
    }
    $clusters[$currentCluster] += ,$currentSequence
}

# Crear el archivo de salida
New-Item -Path $outputFile -ItemType File -Force | Out-Null

# Tomar 50 secuencias al azar de cada cluster y escribirlas en el archivo de salida
foreach ($cluster in $clusters.Keys) {
    $sequences = $clusters[$cluster]
    $randomSequences = $sequences | Get-Random -Count ([Math]::Min(50, $sequences.Count))
    foreach ($sequence in $randomSequences) {
        Add-Content -Path $outputFile -Value $sequence
    }
}

Write-Output "Archivo creado: $outputFile"