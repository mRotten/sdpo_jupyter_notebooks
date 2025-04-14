import os
import comtypes.client


def convert_docx_to_pdf(input_path, output_path=None):
    word = comtypes.client.CreateObject('Word.Application')
    word.Visible = False  # Run Word in the background

    try:
        doc = word.Documents.Open(os.path.abspath(input_path))
        if output_path is None:
            output_path = input_path.replace(".docx", ".pdf")

        doc.SaveAs(os.path.abspath(output_path), FileFormat=17)  # 17 = wdFormatPDF
        doc.Close()
        return output_path
    except Exception as e:
        print(f"Error: {e}")
    finally:
        word.Quit()


if __name__ == "__main__":
    input_docx = "input.docx"  # Change this to your file path
    output_pdf = "output.pdf"  # Change this to desired output path
    convert_docx_to_pdf(input_docx, output_pdf)