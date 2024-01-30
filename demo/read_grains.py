import openpyxl

def read_grains_from_xlsm(filename):
    wb = openpyxl.load_workbook(filename)
    sheet = wb['Sheet1']  # Replace 'Sheet1' with the name of your sheet

    centX = []
    centY = []
    radius = []

    for row in sheet.iter_rows(values_only=True):
        centX.append(row[0])
        centY.append(row[1])
        radius.append(row[2])

    return centX, centY, radius
