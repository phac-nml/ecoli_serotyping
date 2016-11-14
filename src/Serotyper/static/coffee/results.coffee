window.onload = ->

  download = (csv, filename) ->
    csvFile = new Blob([csv], {type: 'text/csv'})
    link = document.createElement('a')
    link.download = filename
    link.href = window.URL.createObjectURL(csvFile)
    link.style.display = 'none'
    document.body.appendChild(link)
    link.click()


  toCSV = (filename) ->
    rows = document.querySelectorAll('table tr')
    file = []

    for r in rows
      row =[]
      cols = r.querySelectorAll('th, td')
      for c in cols
        col = c.innerText
        col = col.split('\n').join('; ')
        row.push(col)
      file.push(row.join(','))

    download(file.join('\r\n'), filename)

  document.querySelector('#download-button').addEventListener('click', (e) ->
      toCSV('SerotypeResults.csv'))